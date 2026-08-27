rule call_cpg:
    group: "mapping"
    input: "output/mapping/{group}_{specimen}.sorted.merged.bam"
    output: "output/cpg_calls/{group}_{specimen}.combined.bed.gz"
    params: prefix = "output/cpg_calls/{group}_{specimen}"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 10
    conda: "../envs/environment.yml"
    shell:
        """
        aligned_bam_to_cpg_scores  --threads 10 --output-prefix {params.prefix} --bam {input}
        """
