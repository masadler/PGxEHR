library(dplyr)
library(tidyr)
library(stringr)

cov_df = read.csv(snakemake@input[["cov"]], sep = '\t')

cohorts = c("statin_LDL", "statin_TC", "statin_HDL", "metformin_HbA1c", 
            "antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
            "antihpt_betablocker_SBP", "betablocker_HR")

colnames <- c("eid", "trait", "filter", "reason")
result_df <- data.frame(matrix(nrow = 0, ncol = length(colnames)))
names(result_df) <- colnames

i = 1
for (cohort in cohorts){
    for (filter in c("s", "l")){
        print(cohort)
        rm_df = read.csv(paste0("output/PGx_phenotype/drug_", cohort, "_", filter, "_1_PGx_removed_individuals.tsv"), sep = '\t')
        rm_df = merge(rm_df, cov_df)
        df = read.csv(paste0("output/PGx_phenotype/drug_", cohort, "_", filter, "_1_PGx_phenotype.tsv"), sep = '\t')

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

        df$reason = "cohort"
        df[df$baseline_measure < l, "reason"] = "low baseline"

        df_cohort = df
        df_cohort = df_cohort[df_cohort$baseline_measure >= l,]
        med_count = df_cohort %>% count(drug_start)
        drugs = as.vector(med_count[med_count$n >= 20, "drug_start"][[1]])

        df[!(df$drug_start %in% drugs) & (df$reason == "cohort"), "reason"] = "outsider drug"

        rm_df = rbind(rm_df[, c("eid", "reason")], df[, c("eid", "reason")])
        rm_df$trait = cohort
        rm_df$filter = filter

        result_df = rbind(result_df, rm_df[, colnames])

        i = i + 1
    }
}

result_df[result_df$reason == "combination medication", "reason"] = "combination therapy"
result_df[result_df$reason == "other antihpt drugs", "reason"] = "other antihpt therapy"
result_df[result_df$reason == "prior antihpt drug", "reason"] = "prior antihpt use"

# aggregate
agg_df <- result_df %>% 
          group_by(trait, filter, reason) %>%
          summarize(count = n())

write.table(agg_df, snakemake@output[[1]], quote = F, row.names = F, sep = '\t')
