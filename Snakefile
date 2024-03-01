import pandas as pd

chroms = [str(i) for i in range(1,23)]

medications = ["statin", "metformin", "betablocker"]
antihpt_medications = ["all", "ACEi", "thiazide_diuretics", "CCB", "betablocker"]
measures = ["LDL", "TC", "HDL", "HbA1c", "SBP", "HR"]

measure_dic = {"LDL": "30780", "TC": "30690", "HDL": "30760", "HbA1c": "30750", "SBP": "4080", "HR": "102"}
measure_file_dic = {"LDL": "ukb672734.csv", "TC": "ukb672734.csv", "HDL": "ukb672734.csv", "HbA1c": "ukb672734.csv",
                    "SBP": "ukb677172.csv", "HR": "ukb677172.csv"}

medication_measure_pairs = ["statin_LDL", "statin_HDL", "statin_TC", "metformin_HbA1c", "betablocker_HR",
                             "antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_thiazide_diuretics_SBP", 
                              "antihpt_betablocker_SBP"]
# "antihpt_CCB_SBP": < 500 samples for stringent filtering

medication_measure_pairs_all = ["statin_LDL", "statin_HDL", "statin_TC", "metformin_HbA1c",
                             "antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_thiazide_diuretics_SBP", "antihpt_CCB_SBP",
                              "antihpt_betablocker_SBP", "betablocker_HR"]

control_medication_measures = ["statin_LDL", "statin_HDL", "statin_TC", "metformin_HbA1c", "betablocker_HR", 
                               "antihpt_all_SBP"] # traits for medication-free individuals

filters = ["s", "l"] # stringent, lenient
n_measures = ["1", "m"] # number of measures, single or multiple
PGx_outcome_phenos = ["pheno_diff", "pheno_log_diff", "pheno_prop", "pheno_diff_adj"] # assessed PGx outcome variations
snp_chrom_dic = {"rs7412": 19, "rs4149056": 12, "rs4970837": 1, "rs11076175": 16, "rs6511720": 19, "rs10455872": 6}

ukbb_path_pheno = "" # Path to UKBB phenotypes
ukbb_path_geno = "" # Path to UKBB genotypes
ukbb_path_org = "" # Path to UKBB organization files
ukbb_path_processed = "" # Path to UKBB processed  files
ref_path = "" # Path to reference panel
gwas_processed = "" # Path to GWAS summary statistics for genetic correlations
gwas_processed2 = "" # Path to GWAS summary statistics for genetic correlations
ldsc_path = "" # Path to ldsc files (LD scores)

ruleorder: prepare_PGx_phenotype_antihpt > prepare_PGx_phenotype
ruleorder: prepare_control_trait_PGx_GWAS > prepare_trait_PGx_GWAS_antihpt > prepare_trait_PGx_GWAS
ruleorder: prepare_control_trait_PGx_INV_GWAS > prepare_trait_PGx_INV_GWAS_antihpt > prepare_trait_PGx_INV_GWAS

rule all:
    input:
        expand("output/medication/{medication}_entries_with_dose.tsv", medication = medications),
        expand("output/measures/{measure}_measures.tsv", measure = measures),
        expand("output/PGx_phenotype/drug_{medication_measure}_{filter}_{n_meas}_PGx_phenotype.tsv", medication_measure = medication_measure_pairs, filter = filters, n_meas = n_measures),
        expand("output/PGx_phenotype/drug_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_phenotype.tsv", medication = antihpt_medications, filter = filters, n_meas = n_measures),
        "output/PGx_phenotype/drug_cohort_statistics.tsv",
        "output/PGx_phenotype/drug_cohort_statistics_filtering.tsv",
        expand("output/GWAS/drug_{medication_measure}_{filter}_{n_meas}_PGx_phenotype_GWAS.txt", medication_measure = medication_measure_pairs, filter = filters, n_meas = n_measures),
        expand("output/GWAS/drug_{medication_measure}_l_m_PGx_phenotype_INV_GWAS.txt", medication_measure = medication_measure_pairs_all),
        expand("output/GWAS/drug_antihpt_CCB_SBP_l_{n_meas}_PGx_phenotype_GWAS.txt", n_meas = n_measures),
        expand("output/GWAS/control_{medication_measure}_l_1_PGx_phenotype_GWAS.txt", medication_measure = control_medication_measures),
        expand("output/GWAS/control_{medication_measure}_l_1_PGx_phenotype_INV_GWAS.txt", medication_measure = control_medication_measures),
        expand("output/GWAS/drug_{medication_measure}_{filter}_{n_meas}_GWAS_{out_pheno}.tsv", medication_measure = medication_measure_pairs, filter = filters, n_meas = n_measures, out_pheno = PGx_outcome_phenos),
        expand("output/GWAS/drug_antihpt_CCB_SBP_l_{n_meas}_GWAS_{out_pheno}.tsv", n_meas = n_measures, out_pheno = PGx_outcome_phenos),
        expand("output/GWAS/control_{medication_measure}_l_1_GWAS_{out_pheno}.tsv", medication_measure = control_medication_measures, out_pheno = PGx_outcome_phenos),
        expand("output/signal/drug_{medication_measure}_{filter}_{n_meas}_GWAS_{out_pheno}_snp_clumped.txt", medication_measure = medication_measure_pairs, filter = filters, n_meas = n_measures, out_pheno = PGx_outcome_phenos),
        expand("output/signal/drug_antihpt_CCB_SBP_l_{n_meas}_GWAS_{out_pheno}_snp_clumped.txt", n_meas = n_measures, out_pheno = PGx_outcome_phenos),
        expand("output/signal/control_{medication_measure}_l_1_GWAS_{out_pheno}_snp_clumped.txt", medication_measure = control_medication_measures, out_pheno = PGx_outcome_phenos),
        "output/Manhattan/drug_statin_LDL_l_m_GWAS_pheno_diff_manhattan.png",
        "output/Manhattan/control_statin_LDL_l_1_GWAS_pheno_diff_manhattan.png",
        expand("output/ldsc/drug_{medication_measure}_l_m_GWASformatted_{out_pheno}.sumstats.gz", medication_measure = medication_measure_pairs_all, out_pheno = ["pheno_diff", "pheno_log_diff"]),
        #"output/ldsc/genetic_correlations.tsv", # GenomicSEM not running inside snakemake as expected, run this step outside of snakemake
        "output/PRS/PRS_correlations.tsv",
        expand("output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baselineres_PRS_absdiffresponse_quintile_stratification.tsv", drug = ["drug"], medication = ["statin"], measure = ["LDL"], filter = ["l"], n_meas = ["m"]),
        expand("output/stratification/{snp}_dosage.tsv", snp = snp_chrom_dic.keys()),
        "output/stratification/stratification_analysis_levels.tsv",
        "output/stratification/stratification_statin_LDL_baseline_PRS_APOE_raw.tsv",
        
rule extract_medication_use:
    input:
        read_v2 = "data/all_lkps_maps_v3_read_v2_drugs.txt",
        bnf_all = "data/all_lkps_maps_v3_bnf_lkp.txt",
        gp_scripts = f'{ukbb_path_pheno}' + "gp_scripts.txt.gz",
        bnf = "data/{medication}_bnf_keyword_search.tsv",
        read2 = "data/{medication}_classification_read2.tsv"
    output:
        no_amount = "output/medication/{medication}_no_amount.tsv",
        no_dose = "output/medication/{medication}_no_dose.tsv",
        dose = "output/medication/{medication}_entries_with_dose.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/extract_primary_care/extract_{wildcards.medication}_use.R"

rule extract_antihypertensive_betablocker_use:
    input:
        read_v2 = "data/all_lkps_maps_v3_read_v2_drugs.txt",
        bnf_all = "data/all_lkps_maps_v3_bnf_lkp.txt",
        gp_scripts = f'{ukbb_path_pheno}' + "gp_scripts.txt.gz",
        bnf = "data/betablocker_bnf_keyword_search.tsv",
        read2_antih = "data/antihypertensive_classification_read2.tsv",
        read2_bb = "data/betablocker_classification_read2.tsv"
    output:
        no_amount = "output/medication/antihypertensive_betablocker_no_amount.tsv",
        no_dose = "output/medication/antihypertensive_betablocker_no_dose.tsv",
        dose = "output/medication/antihypertensive_betablocker_entries_with_dose.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/extract_primary_care/extract_antihypertensive_betablocker_use.R"

rule extract_measures:
    input:
        gp_clinical = f'{ukbb_path_pheno}' + "gp_clinical.txt.gz",
        codes = "data/{measure}_codes.tsv",
        pheno = lambda wildcards: f'{ukbb_path_pheno}' + "{}".format(measure_file_dic[wildcards.measure]),
        inst = f'{ukbb_path_pheno}' + "ukb672662.csv"
    params:
        m = lambda wildcards: measure_dic[wildcards.measure]
    output:
        "output/measures/{measure}_measures.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/extract_primary_care/extract_measures.R"

rule filter_samples:
    input:
        ukb_samples = f'{ukbb_path_pheno}' + "ukb_sqc_v2.txt",
        sample_order = f'{ukbb_path_geno}' +  "plink/ukb1638_cal_chr1_v2_s488366.fam",
        withdrawn = f'{ukbb_path_org}' + "w16389_2023-04-25.csv"
    output:
        "data/UKBB_white_British_samples.txt"
    threads: 1
    resources:
        mem_mb = 2000
    script:
        "scripts/GWAS/filter_samples.R"

rule prepare_covariate_file:
    input:
        cov = f'{ukbb_path_processed}' + "UKBB_covariates.csv", # file with sex, age, PC1-40
        file = f'{ukbb_path_pheno}' + "ukb21067.csv", # phenofile with bmi
        withdrawn = f'{ukbb_path_org}' + "w16389_2023-04-25.csv",
        included = "data/UKBB_white_British_samples.txt"
    params:
        bmi = "21001-0.0",
        year = "34-0.0",
        month = "52-0.0"
    output:
        "output/UKBB_filtered_covariates.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/GWAS/prepare_covariate_file.py"

rule prepare_PGx_phenotype:
    input:
        cov = "output/UKBB_filtered_covariates.tsv",
        medication = "output/medication/{medication}_entries_with_dose.tsv",
        measure = "output/measures/{measure}_measures.tsv",
        gp_clinical = f'{ukbb_path_pheno}' + "gp_clinical.txt.gz",
        gp_scripts = f'{ukbb_path_pheno}' + "gp_scripts.txt.gz"
    threads: 1
    resources:
        mem_mb = 3000
    output:
        "output/PGx_phenotype/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype.tsv",
        "output/PGx_phenotype/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_removed_individuals.tsv"
    script:
        "scripts/prepare_PGx_phenotype/prepare_{wildcards.drug}_{wildcards.medication}_PGx_phenotype.py"

rule prepare_PGx_phenotype_antihpt:
    input:
        cov = "output/UKBB_filtered_covariates.tsv",
        medication = "output/medication/antihypertensive_betablocker_entries_with_dose.tsv",
        measure = "output/measures/SBP_measures.tsv",
        gp_clinical = f'{ukbb_path_pheno}' + "gp_clinical.txt.gz",
        gp_scripts = f'{ukbb_path_pheno}' + "gp_scripts.txt.gz"
    threads: 1
    resources:
        mem_mb = 3000
    output:
        "output/PGx_phenotype/{drug}_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_phenotype.tsv",
        "output/PGx_phenotype/{drug}_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_removed_individuals.tsv"
    script:
        "scripts/prepare_PGx_phenotype/prepare_{wildcards.drug}_antihpt_SBP_PGx_phenotype.py"

rule cohort_statistics:
    input:
        expand("output/PGx_phenotype/drug_{medication_measure}_{filter}_1_PGx_phenotype.tsv", medication_measure = medication_measure_pairs, filter = filters)
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/PGx_phenotype/drug_cohort_statistics.tsv",
        "output/PGx_phenotype/control_cohort_statistics.tsv"
    script:
        "scripts/prepare_PGx_phenotype/cohort_statistics.R"

rule cohort_statistics_filtering:
    input:
        cov = "output/UKBB_filtered_covariates.tsv",
        cohorts = expand("output/PGx_phenotype/drug_{medication_measure}_{filter}_1_PGx_phenotype.tsv", medication_measure = medication_measure_pairs, filter = filters)
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/PGx_phenotype/drug_cohort_statistics_filtering.tsv"
    script:
        "scripts/prepare_PGx_phenotype/cohort_statistics_filtering.R"

rule prepare_trait_PGx_GWAS:
    input:
        pheno = "output/PGx_phenotype/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype.tsv"
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/GWAS/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype_GWAS.txt"
    script:
        "scripts/GWAS/prepare_{wildcards.drug}_{wildcards.medication}_PGx_phenotype_GWAS.R"

rule prepare_trait_PGx_INV_GWAS:
    input:
        pheno = "output/PGx_phenotype/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype.tsv"
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/GWAS/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype_INV_GWAS.txt"
    script:
        "scripts/GWAS/prepare_{wildcards.drug}_{wildcards.medication}_PGx_phenotype_INV_GWAS.R"


rule prepare_trait_PGx_GWAS_antihpt:
    input:
        pheno = "output/PGx_phenotype/{drug}_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_phenotype.tsv",
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/GWAS/{drug}_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_phenotype_GWAS.txt"
    script:
        "scripts/GWAS/prepare_drug_antihpt_SBP_PGx_phenotype_GWAS.R"

rule prepare_trait_PGx_INV_GWAS_antihpt:
    input:
        pheno = "output/PGx_phenotype/{drug}_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_phenotype.tsv",
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/GWAS/{drug}_antihpt_{medication}_SBP_{filter}_{n_meas}_PGx_phenotype_INV_GWAS.txt"
    script:
        "scripts/GWAS/prepare_drug_antihpt_SBP_PGx_phenotype_INV_GWAS.R"

rule prepare_control_trait_PGx_GWAS:
    input:
        pheno = "output/PGx_phenotype/control_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype.tsv"
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/GWAS/control_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype_GWAS.txt"
    script:
        "scripts/GWAS/prepare_control_PGx_phenotype_GWAS.R"

rule prepare_control_trait_PGx_INV_GWAS:
    input:
        pheno = "output/PGx_phenotype/control_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype.tsv"
    threads: 1
    resources:
        mem_mb = 1000
    output:
        "output/GWAS/control_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype_INV_GWAS.txt"
    script:
        "scripts/GWAS/prepare_control_PGx_phenotype_INV_GWAS.R"

rule run_gwas_step1:
    input:
        geno = f'{ukbb_path_processed}' + "geno/_001_ukb_cal_allChrs.bim",
        pheno = "output/GWAS/{trait}_PGx_phenotype_GWAS.txt",
        snps = f'{ukbb_path_processed}' + "geno/qc_pass_MAC100_ICLD_LCR.snplist",
        samples = f'{ukbb_path_processed}' + "geno/qc_pass_MAC100.id"
    threads: 8
    resources:
        mem_mb = 8000
    params:
        geno = f'{ukbb_path_processed}' +  "geno/_001_ukb_cal_allChrs",
        out = "output/GWAS/{trait}_out",
        prefix = "output/GWAS/{trait}" 
    output:
        temp("output/GWAS/{trait}_out.log"),
        temp(expand("output/GWAS/{{trait}}_out_{n}.loco", n = list(range(1, len(PGx_outcome_phenos) + 1)))),
        temp("output/GWAS/{trait}_out_pred.list")  
    shell:
        """
        export LC_ALL=en_US.UTF-8
        export LC_CTYPE=UTF-8
        export LANG=en_US.UTF-8
        regenie_v3.2.6  \
        --step 1 \
        --bed {params.geno} \
        --extract {input.snps} \
        --keep {input.samples} \
        --phenoFile {input.pheno} \
        --bsize 1000 \
        --lowmem \
        --lowmem-prefix {params.prefix} \
        --out {params.out} \
        --threads {threads}
        """

rule run_gwas_step2:
    input:
        step1 = "output/GWAS/{trait}_out_pred.list",
        loco = expand("output/GWAS/{{trait}}_out_{n}.loco", n = list(range(1, len(PGx_outcome_phenos) + 1))),
        pred = "output/GWAS/{trait}_out_pred.list",
        genotypes = f'{ukbb_path_geno}' + "imp/v3/pgen/ukb22828_c{chr}_b0_v3.pgen",
        pheno = "output/GWAS/{trait}_PGx_phenotype_GWAS.txt",
        #sample = f'{ukbb_path_processed}' + "geno/ukb1638_imp_chr1_v2_s487398_sexinfo.sample"
    threads: 8
    resources:
        mem_mb = 4000
    params:
        genotypes = f'{ukbb_path_geno}' + "imp/v3/pgen/ukb22828_c{chr}_b0_v3",
        out = "output/GWAS/UKB_{trait}_chr{chr}"
    output:
        temp(expand("output/GWAS/UKB_{{trait}}_chr{{chr}}_{out_pheno}.regenie", out_pheno = PGx_outcome_phenos)),
        temp("output/GWAS/UKB_{trait}_chr{chr}.log")
    shell:
        """
        export LC_ALL=en_US.UTF-8
        export LC_CTYPE=UTF-8
        export LANG=en_US.UTF-8
        regenie_v3.2.6 \
        --step 2 \
        --pgen {params.genotypes} \
        --ref-first \
        --phenoFile {input.pheno} \
        --chr {wildcards.chr} \
        --approx \
        --pred {input.step1} \
        --pThresh 0.01 \
        --bsize 400 \
        --out {params.out} \
        --threads {threads}
        """

rule aggregate_gwas:
    input:
        files = expand("output/GWAS/UKB_{{trait}}_chr{chr}_{{out_pheno}}.regenie", chr = chroms)
    threads: 1
    resources:
        mem_mb = 15000
    output:
        "output/GWAS/{trait}_GWAS_{out_pheno}.tsv"
    script:
        "scripts/GWAS/aggregate_gwas.py"

rule format_GWAS_regenie:
    input:
        file = "output/GWAS/{trait}_GWAS_{out_pheno}.tsv",
        ref_snps = f'{ref_path}' + "uk10k.autosomal.bim"
    threads: 1
    resources:
        mem_mb = 5000
    output:
        "output/GWAS/{trait}_GWASformatted_{out_pheno}.tsv"
    script:
        "scripts/GWAS/format_GWAS_regenie.py"

rule prune_GWAS:
    input:
        gwas = "output/GWAS/{trait}_GWASformatted_{out_pheno}.tsv",
        bfile = f'{ref_path}' + "uk10k.autosomal.bed"
    params:
        bfile = f'{ref_path}' + "uk10k.autosomal",
        out = "output/signal/{trait}_GWAS_{out_pheno}_snp_clump"
    threads: 1
    resources:
        mem_mb = 5000
    output:
        sig_snps = temp("output/signal/{trait}_GWAS_{out_pheno}_sig_snps.txt"),
        clump_input = temp("output/signal/{trait}_GWAS_{out_pheno}_snp_clump_snp_pval.txt"),
        clump_output = temp("output/signal/{trait}_GWAS_{out_pheno}_snp_clump.clumped"),
        clump_log = temp("output/signal/{trait}_GWAS_{out_pheno}_snp_clump.log"),
        clump_nosex = temp("output/signal/{trait}_GWAS_{out_pheno}_snp_clump.nosex"),
        clumped_snps = "output/signal/{trait}_GWAS_{out_pheno}_snp_clumped.txt"
    shell:
        """
        awk '{{if ($9 <= 5e-8) {{ print $1 }} }}' {input.gwas} > {output.sig_snps}

        if [ -s {output.sig_snps} ]
        then
            awk {{'print $1"\t"$9'}} {input.gwas} > {output.clump_input}

            plink --bfile {params.bfile} --extract {output.sig_snps} --clump {output.clump_input} --clump-r2 0.001 --clump-field 'p' --out {params.out}

            awk '$3 ~ /^rs/ {{print $3 }}' {output.clump_output} > {output.clumped_snps}
        else
            touch {output.clump_input}
            touch {output.clump_output}
            touch {output.clump_log}
            touch {output.clump_nosex}
            touch {output.clumped_snps}
        fi
        """

rule manhattan_plots:
    input:
        gwas = "output/GWAS/drug_{trait}_GWASformatted_{out_pheno}.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    output:
        "output/Manhattan/drug_{trait}_GWAS_{out_pheno}_manhattan.png"
    script:
        "scripts/Figures/manhattan_plots.R"

rule manhattan_plots_control:
    input:
        gwas = "output/GWAS/control_{trait}_GWASformatted_{out_pheno}.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    output:
        "output/Manhattan/control_{trait}_GWAS_{out_pheno}_manhattan.png"
    script:
        "scripts/Figures/manhattan_plots_controls.R"

rule extract_PRS:
    input:
        file = f'{ukbb_path_pheno}' + "ukb676964.csv"
    output:
        "output/PRS/prs_measures.tsv"
    threads: 1
    resources:
        mem_mb = 10000
    script:
        "scripts/PRS/extract_PRS.R"

rule PRS_heatmap:
    input:
        prs = "output/PRS/prs_UKBB_measures.tsv",
        pheno = "output/PGx_phenotype/{drug}_{medication}_{measure}_{filter}_{n_meas}_PGx_phenotype.tsv"
    output:
        absdiff_baseres_heatmap = "output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baselineres_PRS_absdiffresponse_quintile_stratification.tsv",
        absdiff_base_heatmap = "output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baseline_PRS_absdiffresponse_quintile_stratification.tsv",
        ratio_baseres_heatmap = "output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baselineres_PRS_ratioresponse_quintile_stratification.tsv",
        ratio_base_heatmap = "output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baseline_PRS_ratioresponse_quintile_stratification.tsv",
        logdiff_baseres_heatmap = "output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baselineres_PRS_logdiffresponse_quintile_stratification.tsv",
        logdiff_base_heatmap = "output/PRS/{drug}_{medication}_{measure}_{filter}_{n_meas}_baseline_PRS_logdiffresponse_quintile_stratification.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/PRS/PRS_heatmap_{wildcards.drug}_{wildcards.medication}.R"

rule PRS_correlation:
    input:
        prs = "output/PRS/prs_measures.tsv",
        prs_ukbb = "output/PRS/prs_UKBB_measures.tsv",
        phenos = expand("output/PGx_phenotype/drug_{medication_measure}_l_m_PGx_phenotype.tsv", medication_measure = medication_measure_pairs)
    output:
        "output/PRS/PRS_correlations.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/PRS/PRS_correlation.R"

rule munge_sumstats:
    input:
        snps = f'{ldsc_path}' + "w_hm3.snplist",
        base_gwas = expand(f'{gwas_processed}' + "{trait}_gwas_summary_uk10kck.ma", trait = ["CAD_consortium", "HbA1c", "T2D", "SBP_meta", "PR"]),
        glgc_gwas = expand(f'{gwas_processed2}' + "{trait}_summary_stats.tsv", trait = ["LDL_EUR_2021", "HDL_EUR_2021", "TC_EUR_2021"]),
        pgx_gwas = expand("output/GWAS/drug_{medication_measure}_l_m_GWASformatted_{out_pheno}.tsv", medication_measure = medication_measure_pairs_all, out_pheno = ["pheno_diff", "pheno_log_diff"]),
        control_gwas = expand("output/GWAS/control_{medication_measure}_l_1_GWASformatted_{out_pheno}.tsv", medication_measure = control_medication_measures, out_pheno = ["pheno_diff", "pheno_log_diff"])
    output:
        expand("output/ldsc/{trait}.sumstats.gz", trait = ["CAD_consortium", "HbA1c", "T2D", "SBP_meta", "PR"]), # "LDL_EUR_2021", "HDL_EUR_2021", "TC_EUR_2021"
        expand("output/ldsc/drug_{medication_measure}_l_m_GWASformatted_{out_pheno}.sumstats.gz", medication_measure = medication_measure_pairs_all, out_pheno = ["pheno_diff", "pheno_log_diff"]),
        expand("output/ldsc/control_{medication_measure}_l_1_GWASformatted_{out_pheno}.sumstats.gz", medication_measure = control_medication_measures, out_pheno = ["pheno_diff", "pheno_log_diff"])
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/correlation/munge_sumstats.R"

rule genetic_correlations:
    input:
        expand("output/ldsc/{trait}.sumstats.gz", trait = ["CAD_consortium", "HbA1c", "T2D", "SBP_meta", "PR", "LDL_EUR_2021", "HDL_EUR_2021", "TC_EUR_2021"]),
        expand("output/ldsc/drug_{medication_measure}_l_m_GWASformatted_{out_pheno}.sumstats.gz", medication_measure = medication_measure_pairs_all, out_pheno = ["pheno_diff", "pheno_log_diff"]),
        expand("output/ldsc/control_{medication_measure}_l_1_GWASformatted_{out_pheno}.sumstats.gz", medication_measure = control_medication_measures, out_pheno = ["pheno_diff", "pheno_log_diff"])
    params:
        ld_path = f'{ldsc_path}' + "eur_w_ld_chr/"
    output:
        "output/ldsc/genetic_correlations.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    script:
        "scripts/correlation/genetic_correlations.R"

rule extract_stratification:
    input:
        geno = lambda wildcards: f'{ukbb_path_geno}' + "imp/_001_ukb_imp_chr{}_v2.bgen".format(snp_chrom_dic[wildcards.snp]),
        sample_order = f'{ukbb_path_geno}' + "imp/ukb1638_imp_chr1_v2_s487398.sample",
        cov = "output/UKBB_filtered_covariates.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    output:
        "output/stratification/{snp}_dosage.tsv"
    script:
        "scripts/stratification/extract_stratification.R"

rule stratification_analysis:
    input:
        expand("output/stratification/{snp}_dosage.tsv", snp = snp_chrom_dic.keys())
    threads: 1
    resources:
        mem_mb = 5000
    output:
        "output/stratification/stratification_analysis_levels.tsv",
        "output/stratification/stratification_analysis_slopes.tsv"
    script:
        "scripts/stratification/stratification_analysis.R"

rule stratification_analysis_APOE:
    input:
        snp = "output/stratification/rs7412_dosage.tsv",
        pheno = "output/PGx_phenotype/drug_statin_LDL_l_m_PGx_phenotype.tsv",
        prs_ukbb = "output/PRS/prs_UKBB_measures.tsv"
    threads: 1
    resources:
        mem_mb = 5000
    output:
        "output/stratification/stratification_statin_LDL_baseline_PRS_APOE_raw.tsv",
        "output/stratification/stratification_statin_LDL_baseline_PRS_APOE_strata_numerical_values.tsv"
    script:
        "scripts/stratification/stratification_analysis_APOE.R"