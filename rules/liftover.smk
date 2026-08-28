'''
rule bedToFasta:
    input:
        hg38 = config['hg38']['fasta'],
        regions = config['hg38']['regions'],
    output:
        "data/human_regions.fa"
    conda: "../envs/environment.yml"
    shell:
        """
        bedtools getfasta -fi  {input.hg38}  -bed {input.regions} -name  -fo {output}
        """


rule minimap2:
    input:
        reference = config['reference']['fasta'],
        humanRegions = "data/human_regions.fa"
    output:
        "data/hg38_imprinted_regions_to_reference.paf"
    conda: "../envs/environment.yml"
    threads: 16
    shell:
        """
        minimap2 -x asm20 --secondary=no - {input.reference} {input.humanRegions}  {output}
        """
'''

rule liftover: 
    input:
        hg38 = config['hg38']['fasta'],
        regions = config['hg38']['regions'],
        reference = config['reference']['fasta'],
    output:
        liftedRegions = "output/liftover/lifted_imprinted_regions.bed",
        unmappedRegions = "output/liftover/unmapped_imprinted_regions.bed"
    conda: "../envs/environment.yml"
    shell:
        """
        mkdir output/liftover
        if [ "{config[reference][alias]}" != "hg38" ]; then
            python {workflow.basedir}/scripts/minimap2_bed_liftover_preserve_columns.py --bed {input.regions}  --human-fasta {input.hg38} --macaque-ref {input.reference} --output {output.liftedRegions} --unmapped {output.unmappedRegions} --preset asm20  --min-mapq 20
        else
            echo "Skipping liftover because reference alias is hg38"
            cp {input.regions} output/liftover/lifted_imprinted_regions.bed      
        """
