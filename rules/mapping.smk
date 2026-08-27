rule index:
    input:
        refgenome = config['reference']['fasta'],
    output:
        index =  config['reference']['fasta'] + '.mmi'
    conda: "../envs/environment.yml"
    shell:
        """
        pbmm2 index --preset CCS {input.refgenome} {output.index}
        """


rule pbmm2:
    input:
        hifi = "data/{group}/{specimen}/{lane}/{smrtcell}.bam",
        index = config['reference']['fasta'] + '.mmi'
    output:
        temp("output/mapping/temp/{group}_{specimen}/{lane}/{smrtcell}.filt.bam")
    conda: "../envs/environment.yml"
    threads: 16
    shell:
        """
        pbmm2 align {input.index} {input.hifi} {output} -j 16 --preset CCS
        """
