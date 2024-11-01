library(dplyr)
library(tidyr)
library(stringr)

prs_df = read.csv(snakemake@input[["prs"]], sep = '\t')
prs_ukbb_df = read.csv(snakemake@input[["prs_ukbb"]], sep = '\t')
prs_df = merge(prs_df, prs_ukbb_df, by.x = "eid", by.y = "IID")

medication_measure_pairs = c("statin_LDL", "statin_TC", "statin_HDL", "metformin_HbA1c",
                            "antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
                            "antihpt_betablocker_SBP", "betablocker_HR")

colnames = c("trait", "drug", "definition", "PRS", "b", "se", "p", "Unit", "b_std", "se_std", "p_std")
result_df = data.frame(matrix(NA, nrow = 0, ncol = length(colnames)))
names(result_df) = colnames

i = 1
for (pair in medication_measure_pairs){
    print(pair)
    df = read.csv(paste0("output/PGx_phenotype/drug_", pair, "_l_m_PGx_phenotype.tsv"), sep = '\t')

    if (pair == "statin_LDL"){
        drug = "statin"
        unit = "mmol/L"
        l = 2
        covs_filter = "antilipemic_notstatin_prior + drug_start*dose_avg"
        pair_name = "LDL-C response to statin"
        prs_traits = c("LDL_UKBB")
    } else if (pair == "statin_TC"){
        drug = "statin"
        unit = "mmol/L"
        l = 5
        covs_filter = "antilipemic_notstatin_prior + drug_start*dose_avg"
        pair_name = "TC response to statin"
        prs_traits = c("TC_UKBB")
    } else if (pair == "statin_HDL"){
        drug = "statin"
        unit = "mmol/L"
        l = 0.2
        covs_filter = "antilipemic_notstatin_prior + drug_start*dose_avg"
        pair_name = "HDL-C response to statin"
        prs_traits = c("HDL_UKBB")
    } else if (pair == "metformin_HbA1c"){
        drug = "metformin"
        unit = "mmol/mol"
        l = 42
        covs_filter = "BMI*sex + sulfonylureas_prior + sulfonylureas_during + dose_avg"
        pair_name = "HbA1c response to metformin"
        prs_traits = c("HbA1c_UKBB")
    } else if (pair %in% c("antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
        "antihpt_betablocker_SBP")){
        unit = "mmHg"
        l = 120
        if (pair == "antihpt_betablocker_SBP"){
            drug = "beta blocker"
            covs_filter = "drug_start*dose_avg + loop_diuretics_prior"
            pair_name = "SBP response to beta blocker"
        } else {
            drug = "antihypertensives"
            covs_filter = "drug_start*dose_avg + bblocker_prior + loop_diuretics_prior"
            if (pair == "antihpt_all_SBP"){
                pair_name = "SBP response to antihpt (all)"
            } else if (pair == "antihpt_ACEi_SBP"){
                pair_name = "SBP response to antihpt (ACEi)"
            } else if (pair == "antihpt_CCB_SBP"){
                pair_name = "SBP response to antihpt (CCB)"
            } else if (pair == "antihpt_thiazide_diuretics_SBP"){
                pair_name = "SBP response to antihpt (diuretics)"
            }
        }
        prs_traits = c("SBP_UKBB")
    } else if (pair == "betablocker_HR"){
        drug = "beta blocker"
        unit = "beats/min"
        l = 60
        covs_filter = "BMI*sex + drug_start*dose_avg"
        pair_name = "HR response to beta blocker"
        prs_traits = c("HR_UKBB")
    }

    df = df[df$baseline_measure > l,]

    med_count = df %>% count(drug_start)
    drugs = as.vector(med_count[med_count$n >= 20, "drug_start"][[1]])
    df = df[df$drug_start %in% drugs, ]

    if (pair %in% c("antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
        "antihpt_betablocker_SBP")){
        if ((nrow(df[df$loop_diuretics_prior == "yes", ]) < 1)){
            covs_filter = unlist(strsplit(covs_filter, split= ' + loop_diuretics_prior', fixed=TRUE))[1]
        }
    }

    df = merge(df, prs_df)
    df$pheno_diff_out = df$post_measure - df$baseline_measure
    df$pheno_log_diff_out = log(df$post_measure) - log(df$baseline_measure)

    for (prs_trait in prs_traits){
        print(prs_trait)
        df$PRS_std = scale(df[prs_trait])
        for (pheno in c("diff", "log_diff")){

            reg = summary(lm(paste0("pheno_", pheno, "_out ~ PRS_std + sex + age_drug + post_measure_DT + ", covs_filter),
                                data = df, na.action = na.exclude))

            result_df[i, "trait"] = pair_name
            result_df[i, "drug"] = drug
            result_df[i, "PRS"] = prs_trait
            result_df[i, "p"] = reg$coef["PRS_std", "Pr(>|t|)"]
            if (pheno == "diff"){
                result_df[i, "b"] = reg$coef["PRS_std", "Estimate"]
                result_df[i, "se"] = reg$coef["PRS_std", "Std. Error"]
                result_df[i, "definition"] = "post-base"
                result_df[i, "Unit"] = unit
            } else if (pheno == "log_diff"){
                result_df[i, "b"] = reg$coef["PRS_std", "Estimate"]*100
                result_df[i, "se"] = reg$coef["PRS_std", "Std. Error"]*100
                result_df[i, "definition"] = "log(post)-log(base)"
                result_df[i, "Unit"] = "%"
            }

            # scaled regression
            df$pheno_std = scale(df[paste0("pheno_", pheno, "_out")])
            reg = summary(lm(paste0("pheno_std ~ PRS_std + sex + age_drug + post_measure_DT + ", covs_filter),
                                data = df, na.action = na.exclude))
            
            result_df[i, "b_std"] = reg$coef["PRS_std", "Estimate"]
            result_df[i, "se_std"] = reg$coef["PRS_std", "Std. Error"]
            result_df[i, "p_std"] = reg$coef["PRS_std", "Pr(>|t|)"]

            i = i + 1
        }
    }
}

write.table(result_df, snakemake@output[[1]], quote = F, row.names = F, sep = '\t')