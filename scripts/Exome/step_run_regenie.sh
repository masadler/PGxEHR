#!/bin/bash

# regenie - STEP 1

traits=(control_antihpt_all_SBP_l_1 control_betablocker_HR_l_1 control_statin_HDL_l_1 control_statin_LDL_l_1 control_statin_TC_l_1
drug_antihpt_ACEi_SBP_l_m drug_antihpt_CCB_SBP_l_m drug_antihpt_all_SBP_l_m drug_antihpt_betablocker_SBP_l_m drug_antihpt_thiazide_diuretics_SBP_l_m 
drug_betablocker_HR_l_m drug_metformin_HbA1c_l_m drug_statin_HDL_l_m drug_statin_LDL_l_m drug_statin_TC_l_m)
#control_metformin_HbA1c_l_1

for trait in $traits; do
dx upload output/GWAS/${trait}_PGx_phenotype_INV_GWAS.txt --path /PGx_analysis/phenotypeINV/
done

for trait in ${traits}; do
    dx run app-swiss-army-knife \
        -iimage_file="/docker/regenie.tar.gz" \
        -iin="/data/geno/ukb_cal_chr1_22_v2_merged.bed" \
        -iin="/data/geno/ukb_cal_chr1_22_v2_merged.bim" \
        -iin="/data/geno/ukb_cal_chr1_22_v2_merged.fam" \
        -iin="/data/geno/qc_pass_ICLD_LCR_liftover_autosomal.snplist" \
        -iin="/data/geno/qc_pass_temp.id" \
        -iin="/PGx_analysis/phenotypeINV/${trait}_PGx_phenotype_INV_GWAS.txt" \
        -icmd="regenie  \
                --step 1 \
                --nauto 26 \
                --bed ukb_cal_chr1_22_v2_merged \
                --extract qc_pass_ICLD_LCR_liftover_autosomal.snplist \
                --keep qc_pass_temp.id \
                --phenoFile ${trait}_PGx_phenotype_INV_GWAS.txt \
                --bsize 1000 \
                --lowmem \
                --lowmem-prefix ${trait} \
                --out ${trait} \
                --threads 8" \
        --destination "/PGx_analysis/GWAS_step1_INV/" \
        --instance-type mem1_ssd1_v2_x8 --priority normal \
        --name regenie_step1 --tag regenie_step1 \
        -y
done

# regenie - STEP 2 - SKATO

path_to_500kwes_helper_files="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/helper_files/"
genotype_prefix="/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/"

for trait in ${traits}; do
    for CHR in {1..22}; do

        dx run app-swiss-army-knife \
        -iimage_file="/docker/regenie.tar.gz" \
        -iin="/PGx_analysis/GWAS_step1_INV/${trait}_pred.list" \
        -iin="/PGx_analysis/GWAS_step1_INV/${trait}_1.loco" \
        -iin="/PGx_analysis/GWAS_step1_INV/${trait}_2.loco" \
        -iin="/PGx_analysis/GWAS_step1_INV/${trait}_3.loco" \
        -iin="/PGx_analysis/GWAS_step1_INV/${trait}_4.loco" \
        -iin="/data/geno/qc_pass_temp.id" \
        -iin="/PGx_analysis/phenotypeINV/${trait}_PGx_phenotype_INV_GWAS.txt" \
        -iin="${genotype_prefix}ukb23159_c${CHR}_b0_v1.bgen" \
        -iin="${genotype_prefix}ukb23159_c${CHR}_b0_v1.sample" \
        -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.sets.txt.gz" \
        -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
        -iin="/PGx_analysis/LoF_missense_mask.txt" \
        -icmd="regenie  \
                --step 2 \
                --nauto 23 \
                --ref-first \
                --pred ${trait}_pred.list \
                --keep qc_pass_temp.id \
                --sample ukb23159_c${CHR}_b0_v1.sample \
                --bgen ukb23159_c${CHR}_b0_v1.bgen \
                --phenoFile ${trait}_PGx_phenotype_INV_GWAS.txt \
                --set-list ukb23158_500k_OQFE.sets.txt.gz \
                --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
                --mask-def LoF_missense_mask.txt \
                --aaf-bins 0.01 \
                --vc-tests skato \
                --vc-maxAAF 0.01 \
                --bsize 200 \
                --out ${trait}_skato_chr${CHR} \
                --threads 8" \
        --destination "/PGx_analysis/GWAS_step2_INV/" \
        --instance-type mem1_ssd1_v2_x8 --priority normal \
        --name regenie_step2_skato --tag regenie_step2_skato \
        -y

    done
done

# merge files
PGx_outcome_phenos=(pheno_diff pheno_log_diff pheno_prop pheno_diff_adj)

for trait in ${traits}; do
    for out_pheno in ${PGx_outcome_phenos}; do
        merge_cmd="echo $trait ${out_pheno}
        out_file=${trait}_${out_pheno}_INV_skato_WES.tsv
        cp /mnt/project/PGx_analysis/GWAS_step2_INV/${trait}_skato_chr{1..22}_${out_pheno}.regenie .

        echo -e \"CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA\" > \${out_file}

        files=./*.regenie
        for f in \$files
        do
        tail -n+3 \$f | tr \" \" \"\t\" >> \${out_file}
        done
        
        rm *.regenie
        "

        dx run swiss-army-knife \
            -iin="/PGx_analysis/GWAS_step2_INV/${trait}_skato_chr1_${out_pheno}.regenie" \
            -icmd="${merge_cmd}" \
            --tag merge --name merge \
            --instance-type "mem1_ssd1_v2_x2" --priority normal \
            --destination "/PGx_analysis/GWAS_step2_merged_INV/" \
            -y
    done
done

# Download files
for trait in ${traits}; do
    for out_pheno in ${PGx_outcome_phenos}; do
        dx download /PGx_analysis/GWAS_step2_merged_INV/${trait}_${out_pheno}_INV_skato_WES.tsv
    done
done