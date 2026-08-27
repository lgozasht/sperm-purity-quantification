# sperm-purity-quantification
Method to estimate sperm purity from bulk PacBio sperm sequencing data.

## Getting started

To download the pipeline: 
``` https://github.com/lgozasht/sperm-purity-quantification.git ```

You'll  need to generate dataset-specific template_config.yml and template_samples.tsv files before running the pipeline. For each dataset you plan on analyzing with the pipeline, you'll need to generate a new working directory (this directory can be anywhere). Then specify the full path to this directory in the your pipeline config file (e.g. template_config.yml). In this working directory, you also need to make a directory named ```data```. This directory will how the input data for the pipeline with following structure: ```data/{group}/{specimen}/{lane}/{smrtcell}.bam```. These wildcard (specimen, lane, smrtcell) need to be specified in the sample tsv file (e.g. template_samples.tsv). Beware, sample names (or specimens) cannot contain dashes or underscores! The columns in this file provide each of these wildcards:

| specimen | group | lane | smrtcell | 
| :--- | :--- | :--- | :--- |

In the config file you'll need to edit ```workdir```, ```sample_table```, ```reference```, and ```regions```. ```reference``` should just be the full path to the HG38 reference fasta file. ```regions``` should be the path to ```hg38_all_imprinted_supp6_conserved.tsv```.


## Snakemake requirements

```
snakemake-minimal >=8.27
snakemake-executor-plugin-slurm >=0.12.1
snakemake-executor-plugin-slurm-jobstep >=0.2.1
```

## Running the pipeline

Here is an example snakemake command to run the pipeline: 
```snakemake --use-conda --profile PATH/TO/config/snakemake --latency-wait 3 --snakefile /PATH/TO/Snakefile```


## Output

### Summary of reads violating expected methylation signature across imprinted regions for each sample
```{workdir}/output/per_read_cpgs/{group}_{specimen}_final_summaries.tsv```

The last column "proportion_violating" provides the estimated proportion of reads contaminant reads from somatic tissue.

### Per read alignment and methylation summary across imprinted regions 
```{workdir}/output/per_read_cpgs/{group}_{specimen}_per_read_imprinting.tsv  ```




 
