#!/usr/bin/env python3

import argparse
import sys
import pysam
from collections import defaultdict
import statistics

def read_bed(bed_path):
    """
    Read BED file as dict: chrom -> list of (start, end, name)
    BED coordinates are assumed 0-based, half-open.
    """
    regions = defaultdict(list)

    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) >= 4 else f"{chrom}:{start}-{end}"

            regions[chrom].append((start, end, name))

    return regions


def get_modified_c_probs(read, mod_code="m"):
    """
    Return dictionary:
        query_position -> methylation_probability

    Uses pysam's parsed MM/ML tags via read.modified_bases.

    For 5mC, relevant keys are usually:
        ('C', 0, 'm')
        ('C', 1, 'm')

    We combine both strand keys if present.
    Probabilities are converted from 0-255 to 0-1.
    """
    probs = {}

    try:
        mb = read.modified_bases
    except Exception:
        return probs

    if mb is None:
        return probs

    for key, values in mb.items():
        canonical_base, strand, modification = key

        if canonical_base != "C":
            continue

        if modification != mod_code:
            continue

        for qpos, qual in values:
            if qual is None:
                continue
            probs[qpos] = qual / 255.0

    return probs


def is_cpg_site(seq, qpos):
    """
    Return True if qpos is a C in CpG context in the read sequence.

    This uses read sequence context. For reverse-strand alignments, BAM sequence
    handling can be aligner-dependent, but for most PacBio BAMs with MM/ML tags
    this is consistent with query positions in pysam.
    """
    if seq is None:
        return False

    if qpos < 0 or qpos >= len(seq):
        return False

    base = seq[qpos].upper()

    if base != "C":
        return False

    if qpos + 1 >= len(seq):
        return False

    return seq[qpos + 1].upper() == "G"


def is_c_site(seq, qpos):
    if seq is None:
        return False
    if qpos < 0 or qpos >= len(seq):
        return False
    return seq[qpos].upper() == "C"


def process_region(
    bam,
    chrom,
    region_start,
    region_end,
    min_mapq,
    meth_threshold,
    context,
    mod_code,
):
    """
    Process all reads overlapping one BED region.
    """

    for read in bam.fetch(chrom, region_start, region_end):
        if read.is_unmapped:
            continue

        if read.is_secondary or read.is_supplementary:
            continue

        if read.mapping_quality < min_mapq:
            continue

        aln_start = read.reference_start
        aln_end = read.reference_end

        if aln_start is None or aln_end is None:
            continue

        overlap_start = max(region_start, aln_start)
        overlap_end = min(region_end, aln_end)

        if overlap_start >= overlap_end:
            continue

        seq = read.query_sequence
        if seq is None:
            continue

        meth_probs = get_modified_c_probs(read, mod_code=mod_code)

        qpos_in_overlap = []
        site_probs = []

        # aligned_pairs gives query/reference coordinate mapping
        # matches_only=False keeps deletions/skips as None; we ignore those.
        for qpos, rpos in read.get_aligned_pairs(matches_only=False):
            if qpos is None or rpos is None:
                continue

            if rpos < overlap_start or rpos >= overlap_end:
                continue

            if context == "CpG":
                is_site = is_cpg_site(seq, qpos)
            elif context == "C":
                is_site = is_c_site(seq, qpos)
            else:
                raise ValueError("context must be 'CpG' or 'C'")

            if not is_site:
                continue

            p = meth_probs.get(qpos)
            if p is None:
                continue
            qpos_in_overlap.append(qpos)

            # If a cytosine is not listed in MM/ML, treat as unmethylated.
            # This is appropriate for many PacBio 5mC BAMs.
            
            site_probs.append(p)

        if len(qpos_in_overlap) == 0:
            continue

        read_region_start = min(qpos_in_overlap)
        read_region_end = max(qpos_in_overlap) + 1

        n_sites = len(site_probs)
        n_methylated = sum(1 for p in site_probs if p >= meth_threshold)
        #mean_methylation = sum(site_probs) / n_sites
        #mean_methylation = statistics.median(site_probs)
        mean_methylation = n_methylated/n_sites

        yield [
            read.query_name,
            chrom,
            overlap_start,
            overlap_end,
            read_region_start,
            read_region_end,
            n_sites,
            n_methylated,
            mean_methylation,
        ]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-read mean methylation across BED regions "
            "from a PacBio BAM with MM/ML methylation tags."
        )
    )

    parser.add_argument(
        "-b",
        "--bed",
        required=True,
        help="BED file of regions, 0-based half-open.",
    )

    parser.add_argument(
        "-a",
        "--bam",
        required=True,
        help="Indexed BAM/CRAM file containing PacBio methylation tags.",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file.",
    )

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality. Default: 0.",
    )

    parser.add_argument(
        "--meth-threshold",
        type=float,
        default=0.5,
        help=(
            "Probability threshold for counting a site as methylated. "
            "Default: 0.5."
        ),
    )

    parser.add_argument(
        "--context",
        choices=["CpG", "C"],
        default="CpG",
        help=(
            "Which cytosines to count as methylation sites. "
            "Use 'CpG' for CpG cytosines, or 'C' for all cytosines. "
            "Default: CpG."
        ),
    )

    parser.add_argument(
        "--mod-code",
        default="m",
        help=(
            "Modification code to use from MM tags. "
            "For 5mC this is usually 'm'. Default: m."
        ),
    )

    args = parser.parse_args()

    regions = read_bed(args.bed)

    bam = pysam.AlignmentFile(args.bam, "rb")

    with open(args.output, "w") as out:
        header = [
            "read",
            "aligned_chromosome",
            "start_of_alignment_in_region",
            "end_of_alignment_in_region",
            "start_alignment_in_read",
            "end_alignment_in_read",
            "N_sites",
            "N_methylated_sites",
            "mean_methylation_across_alignment",
        ]
        out.write("\t".join(header) + "\n")

        for chrom in regions:
            if chrom not in bam.references:
                print(
                    f"Warning: chromosome {chrom} not found in BAM; skipping.",
                    file=sys.stderr,
                )
                continue

            for region_start, region_end, region_name in regions[chrom]:
                for row in process_region(
                    bam=bam,
                    chrom=chrom,
                    region_start=region_start,
                    region_end=region_end,
                    min_mapq=args.min_mapq,
                    meth_threshold=args.meth_threshold,
                    context=args.context,
                    mod_code=args.mod_code,
                ):
                    out.write("\t".join(map(str, row)) + "\n")

    bam.close()


if __name__ == "__main__":
    main()
