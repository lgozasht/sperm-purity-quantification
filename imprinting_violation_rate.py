#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def read_regions(region_file):
    """
    Read imprinting region file.

    Expected columns:
        chrom start end strand parent gene

    Example:
        chr10 119726042 119829147 + Paternal INPP5F
    """
    regions = defaultdict(list)

    with open(region_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            if len(fields) < 6:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[3]
            parent = fields[4]
            gene = fields[5]

            regions[chrom].append(
                {
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "parent": parent,
                    "gene": gene,
                }
            )

    return regions


def overlaps(a_start, a_end, b_start, b_end):
    """
    Half-open interval overlap:
        [a_start, a_end) overlaps [b_start, b_end)
    """
    return a_start < b_end and b_start < a_end


def classify_read(mean_methylation, threshold):
    if mean_methylation >= threshold:
        return "methylated"
    else:
        return "unmethylated"


def expected_state(parent):
    parent = parent.lower()

    if parent == "maternal":
        return "unmethylated"
    elif parent == "paternal":
        return "methylated"
    else:
        return "unknown"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Report proportion of reads that violate expected methylation "
            "state for paternal/maternal imprinted regions."
        )
    )

    parser.add_argument(
        "-m",
        "--methylation",
        required=True,
        help="Per-read methylation TSV output file.",
    )

    parser.add_argument(
        "-r",
        "--regions",
        required=True,
        help="Imprinted region annotation file.",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file with per-read annotations.",
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Output summary TSV file.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Methylation threshold for classifying reads. Default: 0.5.",
    )

    args = parser.parse_args()

    regions = read_regions(args.regions)

    total_reads = 0
    violating_reads = 0
    matching_reads = 0
    unannotated_reads = 0

    per_gene_stats = defaultdict(lambda: {
        "total": 0,
        "violating": 0,
        "matching": 0,
    })

    per_parent_stats = defaultdict(lambda: {
        "total": 0,
        "violating": 0,
        "matching": 0,
    })

    with open(args.methylation) as infile, open(args.output, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")

        output_fields = reader.fieldnames + [
            "gene",
            "parent",
            "expected_state",
            "observed_state",
            "violates_expectation",
        ]

        writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=output_fields)
        writer.writeheader()

        for row in reader:
            chrom = row["aligned_chromosome"]
            start = int(row["start_of_alignment_in_region"])
            end = int(row["end_of_alignment_in_region"])
            mean_methylation = float(row["mean_methylation_across_alignment"])

            observed = classify_read(mean_methylation, args.threshold)

            overlapping_regions = []

            for region in regions.get(chrom, []):
                if overlaps(start, end, region["start"], region["end"]):
                    overlapping_regions.append(region)

            if len(overlapping_regions) == 0:
                unannotated_reads += 1
                outrow = dict(row)
                outrow.update(
                    {
                        "gene": "NA",
                        "parent": "NA",
                        "expected_state": "NA",
                        "observed_state": observed,
                        "violates_expectation": "NA",
                    }
                )
                writer.writerow(outrow)
                continue

            # If a read overlaps multiple annotated regions, report one row per overlap.
            for region in overlapping_regions:
                parent = region["parent"]
                gene = region["gene"]
                expected = expected_state(parent)

                if expected == "unknown":
                    violates = "NA"
                else:
                    violates = observed != expected

                outrow = dict(row)
                outrow.update(
                    {
                        "gene": gene,
                        "parent": parent,
                        "expected_state": expected,
                        "observed_state": observed,
                        "violates_expectation": violates,
                    }
                )
                writer.writerow(outrow)

                if violates == "NA":
                    continue

                total_reads += 1

                per_gene_stats[gene]["total"] += 1
                per_parent_stats[parent]["total"] += 1

                if violates:
                    violating_reads += 1
                    per_gene_stats[gene]["violating"] += 1
                    per_parent_stats[parent]["violating"] += 1
                else:
                    matching_reads += 1
                    per_gene_stats[gene]["matching"] += 1
                    per_parent_stats[parent]["matching"] += 1

    with open(args.summary, "w") as out:
        out.write("category\tgroup\ttotal_reads\tmatching_reads\tviolating_reads\tproportion_violating\n")

        if total_reads > 0:
            prop_violating = violating_reads / total_reads
        else:
            prop_violating = "NA"

        out.write(
            f"overall\tall\t{total_reads}\t{matching_reads}\t{violating_reads}\t{prop_violating}\n"
        )

        for parent, stats in sorted(per_parent_stats.items()):
            total = stats["total"]
            matching = stats["matching"]
            violating = stats["violating"]
            prop = violating / total if total > 0 else "NA"

            out.write(
                f"parent\t{parent}\t{total}\t{matching}\t{violating}\t{prop}\n"
            )

        for gene, stats in sorted(per_gene_stats.items()):
            total = stats["total"]
            matching = stats["matching"]
            violating = stats["violating"]
            prop = violating / total if total > 0 else "NA"

            out.write(
                f"gene\t{gene}\t{total}\t{matching}\t{violating}\t{prop}\n"
            )

        out.write(
            f"unannotated\tNA\t{unannotated_reads}\tNA\tNA\tNA\n"
        )


if __name__ == "__main__":
    main()
