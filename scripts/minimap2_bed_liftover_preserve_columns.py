#!/usr/bin/env python3

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def read_bed_with_ids(bed_path):
    """
    Read BED and assign each row a unique temporary ID.
    Replace column 4 with the temporary ID.

    Returns:
        rows_by_id: dict temp_id -> original columns with col4 replaced by temp_id
        temp_bed_lines: BED lines with temp_id in col4
    """
    rows_by_id = {}
    temp_bed_lines = []

    with open(bed_path) as f:
        for i, line in enumerate(f, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                raise ValueError(f"BED line has fewer than 3 columns: {line}")

            chrom = fields[0]
            start = fields[1]
            end = fields[2]

            temp_id = f"region_{i}"

            # Ensure at least 4 columns
            if len(fields) < 4:
                fields = fields + [temp_id]
            else:
                fields[3] = temp_id

            rows_by_id[temp_id] = fields
            temp_bed_lines.append("\t".join([chrom, start, end, temp_id]))

    return rows_by_id, temp_bed_lines


def write_temp_bed(temp_bed_lines, out_path):
    with open(out_path, "w") as out:
        for line in temp_bed_lines:
            out.write(line + "\n")


def run_bedtools_getfasta(human_fasta, temp_bed, output_fasta):
    cmd = [
        "bedtools", "getfasta",
        "-fi", human_fasta,
        "-bed", temp_bed,
        "-name",
        "-fo", output_fasta
    ]

    subprocess.run(cmd, check=True)


def run_minimap2_paf(macaque_fasta_or_index, query_fasta, output_paf, preset="asm20"):
    cmd = [
        "minimap2",
        "-x", preset,
        "--secondary=no",
        macaque_fasta_or_index,
        query_fasta
    ]

    with open(output_paf, "w") as out:
        subprocess.run(cmd, check=True, stdout=out)


def clean_query_name(qname):
    """
    bedtools getfasta -name may produce names like:
        region_1::chr1:100-200
    or sometimes just:
        region_1

    This keeps only the temporary region ID.
    """
    return qname.split("::")[0]


def parse_best_paf_hits(paf_path, min_mapq=0):
    """
    Parse PAF and keep one best hit per query.

    Ranking:
        1. highest MAPQ
        2. longest alignment block length
        3. highest number of matching bases

    PAF columns:
        1 query name
        2 query length
        3 query start
        4 query end
        5 strand
        6 target name
        7 target length
        8 target start
        9 target end
        10 residue matches
        11 alignment block length
        12 mapping quality
    """
    best = {}

    with open(paf_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 12:
                continue

            qname = clean_query_name(fields[0])
            strand = fields[4]
            target = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            aln_block_len = int(fields[10])
            mapq = int(fields[11])

            if mapq < min_mapq:
                continue

            hit = {
                "target": target,
                "tstart": tstart,
                "tend": tend,
                "strand": strand,
                "matches": matches,
                "aln_block_len": aln_block_len,
                "mapq": mapq,
            }

            if qname not in best:
                best[qname] = hit
            else:
                old = best[qname]
                old_rank = (old["mapq"], old["aln_block_len"], old["matches"])
                new_rank = (hit["mapq"], hit["aln_block_len"], hit["matches"])

                if new_rank > old_rank:
                    best[qname] = hit

    return best


def write_lifted_bed(rows_by_id, best_hits, output_bed, unmapped_bed=None):
    with open(output_bed, "w") as out:
        unmapped_out = open(unmapped_bed, "w") if unmapped_bed else None

        for temp_id, original_fields in rows_by_id.items():
            if temp_id not in best_hits:
                if unmapped_out:
                    unmapped_out.write("\t".join(original_fields) + "\n")
                continue

            hit = best_hits[temp_id]

            lifted_fields = original_fields.copy()

            # Replace coordinate columns with macaque coordinates
            lifted_fields[0] = hit["target"]
            lifted_fields[1] = str(hit["tstart"])
            lifted_fields[2] = str(hit["tend"])

            # Column 4 can be whatever; use temp_id
            lifted_fields[3] = temp_id

            # Preserve all columns from column 5 onward exactly
            out.write("\t".join(lifted_fields) + "\n")

        if unmapped_out:
            unmapped_out.close()


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Map human BED regions to macaque with minimap2 and preserve the original BED column structure. "
            "Columns 1-3 are replaced with macaque coordinates. Column 4 is replaced with a temporary region ID. "
            "Columns 5 onward are preserved."
        )
    )

    parser.add_argument(
        "--bed",
        required=True,
        help="Input human BED file."
    )

    parser.add_argument(
        "--human-fasta",
        required=True,
        help="Human reference FASTA used to extract input BED sequences."
    )

    parser.add_argument(
        "--macaque-ref",
        required=True,
        help="Macaque FASTA or minimap2 index .mmi."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output lifted macaque BED."
    )

    parser.add_argument(
        "--unmapped",
        default=None,
        help="Optional BED file for unmapped input regions."
    )

    parser.add_argument(
        "--preset",
        default="asm20",
        help="minimap2 preset. Default: asm20."
    )

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to keep. Default: 0."
    )

    args = parser.parse_args()

    rows_by_id, temp_bed_lines = read_bed_with_ids(args.bed)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        temp_bed = tmpdir / "regions.with_temp_ids.bed"
        temp_fasta = tmpdir / "regions.fa"
        temp_paf = tmpdir / "regions_to_macaque.paf"

        write_temp_bed(temp_bed_lines, temp_bed)

        run_bedtools_getfasta(
            human_fasta=args.human_fasta,
            temp_bed=str(temp_bed),
            output_fasta=str(temp_fasta)
        )

        run_minimap2_paf(
            macaque_fasta_or_index=args.macaque_ref,
            query_fasta=str(temp_fasta),
            output_paf=str(temp_paf),
            preset=args.preset
        )

        best_hits = parse_best_paf_hits(
            paf_path=str(temp_paf),
            min_mapq=args.min_mapq
        )

        write_lifted_bed(
            rows_by_id=rows_by_id,
            best_hits=best_hits,
            output_bed=args.output,
            unmapped_bed=args.unmapped
        )


if __name__ == "__main__":
    main()
