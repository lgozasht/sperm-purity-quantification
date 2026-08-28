rule perReadStats:
    input: 
        bamFile = "output/mapping/{group}_{specimen}.sorted.merged.bam",
        bamIndex = "output/mapping/{group}_{specimen}.sorted.merged.bam.bai"
    output: "output/per_read_cpgs/{group}_{specimen}.tsv"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    params:
        regions = "output/liftover/lifted_imprinted_regions.bed"
    threads: 10
    conda: "../envs/environment.yml"
    shell:
        """
        python {workflow.basedir}/scripts/mean_read_methylation_in_bed.py --bed {params.regions}  --bam {input.bamFile}   --output {output}  --mod-code m --context CpG
        """


rule summarizeImprintingViolations:
    input:
        "output/per_read_cpgs/{group}_{specimen}.tsv"
    output:
        perReadSummary = "output/per_read_cpgs/{group}_{specimen}_per_read_imprinting.tsv",
        finalSummary = "output/per_read_cpgs/{group}_{specimen}_final_summaries.tsv"
    params:
        regions =  "output/liftover/lifted_imprinted_regions.bed"
    shell:
        """
        python {workflow.basedir}/scripts/imprinting_violation_rate.py  -m {input} -r {params.regions}  -o {output.perReadSummary} --summary {output.finalSummary} --threshold 0.5
        """
