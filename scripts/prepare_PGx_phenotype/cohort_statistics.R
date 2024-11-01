library(dplyr)
library(tidyr)
library(stringr)

#### drug cohorts
cohorts = c("statin_LDL", "statin_TC", "statin_HDL", "metformin_HbA1c", 
            "antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
            "antihpt_betablocker_SBP", "betablocker_HR")

colnames = c("cohort", "N", "Female", "Age", "BMI", "baseline_time", "baseline", "post_time", "post_level",
             "completeness", "drugs")
result_df = data.frame(matrix(NA, length(cohorts), length(colnames),))
names(result_df) = colnames

i = 1
for (cohort in cohorts){
    for (filter in c("s", "l")){
        print(filter)
        print(cohort)
        df <- read.csv(paste0("output/PGx_phenotype/drug_", cohort, "_", filter, "_1_PGx_phenotype.tsv"), sep = '\t')

        if (cohort == "statin_LDL"){
            unit = "mmol/L"
            l = 2
        } else if (cohort == "statin_TC"){
            unit = "mmol/L"
            l = 5
        } else if (cohort == "statin_HDL"){
            unit = "mmol/L"
            l = 0.2
         } else if (cohort == "metformin_HbA1c"){
            unit = "mmol/mol"
            l = 42
        } else if (cohort %in% c("antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
            "antihpt_betablocker_SBP")){
            unit = "mmHg"
            l = 120
        } else if (cohort == "betablocker_HR"){
            unit = "beats/min"
            l = 60
        }

        df = df[df$baseline_measure >= l,]

        med_count = df %>% count(drug_start)
        drugs = as.vector(med_count[med_count$n >= 20, "drug_start"][[1]])
        print(drugs)
        df = df[df$drug_start %in% drugs, ]

        result_df[i, "cohort"] = cohort
        n = nrow(df)
        result_df[i, "N"] = n
        n_fem = nrow(df[df$sex == 0, ])
        result_df[i, "Female"] = paste0(n_fem, " (", signif(n_fem/n*100, 3), ")")
        result_df[i, "Age"] = paste0(signif(mean(df$age_drug), 3), " (", signif(sd(df$age_drug), 3), ")")
        result_df[i, "BMI"] = paste0(signif(mean(df$BMI, na.rm=TRUE), 3), " (", signif(sd(df$BMI, na.rm=TRUE), 3), ")")
        result_df[i, "baseline_time"] = paste0(signif(mean(df$baseline_DT), 3), " (", signif(sd(df$baseline_DT), 3), ") / ",
        signif(median(df$baseline_DT), 3), "[", signif(quantile(df$baseline_DT, 0.25), 3), ", ", signif(quantile(df$baseline_DT, 0.75), 3),"]")
        result_df[i, "baseline"] = paste0(signif(mean(df$baseline_measure), 3), " (", signif(sd(df$baseline_measure), 3), ") ", unit)
        result_df[i, "post_time"] = paste0(signif(mean(df$post_measure_DT), 3), " (", signif(sd(df$post_measure_DT), 3), ") / ",
        signif(median(df$post_measure_DT), 3), "[", signif(quantile(df$post_measure_DT, 0.25), 3), ", ", signif(quantile(df$post_measure_DT, 0.75), 3), "]")
        result_df[i, "post_level"] = paste0(signif(mean(df$post_measure), 3), " (", signif(sd(df$post_measure), 3), ") ", unit)
        result_df[i, "completeness"] = paste0(signif(mean(df$completeness*100), 3), " (", signif(sd(df$completeness*100), 3), ")")

        drug_stats = ""

        for (d in drugs){
            n_drug = nrow(df[df$drug_start == d, ])
            drug_stats = paste0(drug_stats, d, ": ", n_drug, " (", signif(n_drug/n*100, 3), "); ")
        }

        drug_stats = str_sub(drug_stats, start = 1, end = -3) # remove last semicolon

        print(drug_stats)

        result_df[i, "drugs"] = drug_stats

        i = i + 1
    }
}

write.table(result_df, snakemake@output[[1]], quote = F, row.names = F, sep = '\t')

#### control cohorts

cohorts = c("statin_LDL", "statin_TC", "statin_HDL", "metformin_HbA1c", 
            "antihpt_all_SBP", "betablocker_HR")

colnames = c("cohort", "N", "Female", "Age", "BMI", "measure1", "measure2", "baseline_post_time")
result_df = data.frame(matrix(NA, length(cohorts), length(colnames),))
names(result_df) = colnames

i = 1
for (cohort in cohorts){
    df <- read.csv(paste0("output/PGx_phenotype/control_", cohort, "_l_1_PGx_phenotype.tsv"), sep = '\t')

    if (cohort == "statin_LDL"){
        unit = "mmol/L"
    } else if (cohort == "statin_TC"){
        unit = "mmol/L"
    } else if (cohort == "statin_HDL"){
        unit = "mmol/L"
        } else if (cohort == "metformin_HbA1c"){
        unit = "mmol/mol"
    } else if (cohort == "antihpt_all_SBP"){
        unit = "mmHg"
    } else if (cohort == "betablocker_HR"){
        unit = "beats/min"
    }

    result_df[i, "cohort"] = cohort
    n = nrow(df)
    result_df[i, "N"] = n
    n_fem = nrow(df[df$sex == 0, ])
    result_df[i, "Female"] = paste0(n_fem, " (", signif(n_fem/n*100, 3), ")")
    result_df[i, "Age"] = paste0(signif(mean(df$age_m1), 3), " (", signif(sd(df$age_m1), 3), ")")
    result_df[i, "BMI"] = paste0(signif(mean(df$BMI, na.rm=TRUE), 3), " (", signif(sd(df$BMI, na.rm=TRUE), 3), ")")
    result_df[i, "measure1"] = paste0(signif(mean(df$measure1), 3), " (", signif(sd(df$measure1), 3), ") ", unit)
    result_df[i, "measure2"] = paste0(signif(mean(df$measure2), 3), " (", signif(sd(df$measure2), 3), ") ", unit)
    result_df[i, "baseline_post_time"] = paste0(signif(mean(df$measure_DT), 3), " (", signif(sd(df$measure_DT), 3), ")")

    i = i + 1
}

write.table(result_df, snakemake@output[[2]], quote = F, row.names = F, sep = '\t')