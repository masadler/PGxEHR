#!/bin/bash

plink2 \
  --write-samples \
  --write-snplist \
  --autosome \
  --bfile _001_ukb_cal_allChrs \
  --mind 0.1 \
  --geno 0.01 --maf 0.01 --hwe 1e-15 \
  --exclude range /data/FAC/FBM/DBC/zkutalik/default_sensitive/msadler/snakemake/UKBB_PGx/data/high-LD-regions-hg19-GRCh37.bed \
  --indep-pairwise 1000 100 0.9 \
  --out qc_pass_MAC100_temp

# High-LD regions: https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.bed
awk 'NR == FNR {a[$1]; next} !($1 in a)' /data/FAC/FBM/DBC/zkutalik/default_sensitive/msadler/snakemake/UKBB_PGx/data/ICLD_snps.txt qc_pass_MAC100_temp.prune.in >  qc_pass_MAC100_ICLD_LCR.snplist 