# Inject genotypes into a VCF file

This Snakemake workflow is designed for extract phased genotype information from A VCF file and "inject" into another VCF file. 
It is intended for be used as preparation step for [MSMC2](https://github.com/stschiff/msmc2) that requires input VCF in a specific format, but also requires phasing information. 

The pipeline also calculates statistics from the log files. The statistics are written to `results/replacment_statistics.tsv'.

## Overview

The workflow consists of three main steps:

1. **Download Files (rule `download_files`):**
   - Downloads phased and unphased VCF files for each sample and chromosome combination from a remote server. 
   - This can be deactivated in the config fie. 
   - Output:
     - Unphased VCF file: `test/vcf-unphased/{sample}.{chr}.vcf.gz`
     - Phased VCF file: `test/vcf-phased/{sample}_{chr}.vcf.gz`

2. **Merge Genotypes (rule `merge_genotypes`):**
   - Merges the phased and unphased VCF files for each sample and chromosome.
   - Output:
     - Merged VCF file: `test/vcf-merged/{sample}.{chr}.vcf.gz`
   - Log file: `logs/merge_genotypes/{sample}_{chr}.log`

3. **Get Statistics (rule `get_stats`):**
   - Aggregates statistics from the merge operation for each sample and chromosome.
   - Output:
     - Statistics Overview: `results/overview.tsv`

## Usage

1. Clone the repository:

    ```bash
    git clone https://github.com/mobilegenome/msmc2-prepare.git
    ```

2. Install Snakemake:

    ```bash
    conda install -c bioconda -c conda-forge snakemake
    ```

3. Customize the `config.yaml` file if necessary.

4. Execute the workflow:

    ```bash
    snakemake --cores <num_cores>
    ```

Replace `<num_cores>` with the desired number of CPU cores.

## Configuration

The workflow is configurable via the `config.yaml` file. Adjust the parameters according to your specific setup and requirements.


## Dependencies

- Python 3.x
- Snakemake
- R (for statistics aggregation)

## Scripts

- `workflow/scripts/merge_genotypes.py`: Python script for merging phased and unphased genotypes.
- `workflow/scripts/aggregate_stats.R`: R script for aggregating statistics from the merged genotypes.
