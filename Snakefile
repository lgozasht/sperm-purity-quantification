from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd

configfile: "/global/scratch/users/landen_gozashti/projects/Sperm_diversity/methylation_purity/human/config/pipelineConfigs/human_config.yml"
#configfile: "/global/scratch/users/landen_gozashti/projects/Sperm_diversity/methylation_purity/macaque/config/pipelineConfigs/macaque_config.yml"


workdir: config['workdir'] #Wrong syntax?


### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen']
groups =sample_table["group"]
ref_fasta = config['reference']['fasta']

print(specimens)
print(groups)



# include helper functions
include: "rules/samtools_utils.smk"

include: "rules/mapping.smk"

include: "rules/call_cpgs.smk"

include: "rules/perReadStats.smk"



rule all:
    input:
        # self-alignment: assembly + QC, variant calls through qc_all stage
        expand("output/mapping/{group}_{specimen}.sorted.merged.bam.bai",zip, group = groups, specimen = specimens),
        expand("output/per_read_cpgs/{group}_{specimen}.tsv",zip, group = groups, specimen = specimens),
        expand("output/per_read_cpgs/{group}_{specimen}_final_summaries.tsv",zip, group = groups, specimen = specimens),
        expand("output/cpg_calls/{group}_{specimen}.combined.bed.gz",zip, group = groups, specimen = specimens)
