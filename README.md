# Pharmacogenetic analysis of EHRs

Welcome to the Github accompanying the manuscript **Leveraging large-scale biobank EHRs to enhance pharmacogenetics of cardiometabolic disease medications** which analyses the genetics of drug response for biomarker-medication pairs.

# Usage

This Github contains the workflow pipeline that was used to extract medication records and biomarker measures from electronic health records (EHRs), define pharmacogenetic (PGx) phenotypes and run PGx-GWAS analyses. Code to run whole exome burden tests on the [UK Biobank DNAnexus Research Analysis Platform](https://ukbiobank.dnanexus.com/) is also provided.

The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. This allows for a reproducible analysis and parallel computations - this is especially useful for parallel GWAS compuations. However, individual scripts can also be run without the snakemake workflow manager by replacing snakemake variables by hard-coded input and output paths.

## Input data

Individual-level data from the UK Biobank is needed with access to the primary care records from the following files:

- GP prescription records (datafield \#42039, `gp_scripts.txt`)
- GP clinical event records (datafield \#42040, `gp_clinical.txt`)
- From the `all_lkps_maps_v3.xlsx` downloadable from [Resource 592](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592), the sheet `read_v2_drugs_lkp` called herein: `all_lkps_maps_v3_read_v2_drugs.txt`
- From the `all_lkps_maps_v3.xlsx` downloadable from [Resource 592](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592), the sheet `bnf_lkp` called herein: `all_lkps_maps_v3_bnf_lkp.txt`

## Drug response phenotypes

Drug response phenotypes were derived by:

1) extracting drug prescriptions of the broad medication class and medication of interest (e.g. lipid-lowering medications and statins)

2) extracting biomarker measures from the EHRs and combining them with biomarker measures from the UKBB assessment visits

3) combining these two datasets to extract baseline and post-treatment measures using temporal medication and biomarker information. In addition to the presence of baseline and post-treatment measures in the defined time window relative to medication start, individuals had to pass several other QC steps to be included in a PGx cohort such as consistent drug prescriptions with no change in the drug regimen. 

## Genetic analyses

Genome-wide association analyses (GWAS) were run using regenie in a [2-step procedure](https://rgcgithub.github.io/regenie/). The first step includes a whole genome regression performed on genotyped SNPs passing several QC filters (`scripts/GWAS/QC_SNPs.sh`). The second step includes the testing of genotyped/imputed SNPs in a LOCO scheme.

Rare variant burden tests on whole exome sequencing data was performed on the UK Biobank DNAnexus research analysis platform. Corresponding scripts are in the `scripts/Exome` folder.

## Software requirements

This workflow has been tested with snakemake v7.30.1, python v3.9.13 and R v4.3.2.

# License

MIT License

