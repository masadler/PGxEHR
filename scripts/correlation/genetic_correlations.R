library(GenomicSEM)

setwd("output/ldsc/")

#### genetic correlations ####

ld = snakemake@params[["ld_path"]]

# trait name pairs to test:
traits_1 = c("drug_statin_LDL_l_m_GWASformatted_pheno_diff", "drug_statin_LDL_l_m_GWASformatted_pheno_log_diff",
            "drug_statin_LDL_l_m_GWASformatted_pheno_diff", "drug_statin_LDL_l_m_GWASformatted_pheno_log_diff",
             "LDL_EUR_2021", "LDL_EUR_2021",
            "drug_statin_HDL_l_m_GWASformatted_pheno_diff", "drug_statin_HDL_l_m_GWASformatted_pheno_log_diff",
            "drug_statin_HDL_l_m_GWASformatted_pheno_diff", "drug_statin_HDL_l_m_GWASformatted_pheno_log_diff",
             "HDL_EUR_2021", "HDL_EUR_2021",
            "drug_statin_TC_l_m_GWASformatted_pheno_diff", "drug_statin_TC_l_m_GWASformatted_pheno_log_diff",
            "drug_statin_TC_l_m_GWASformatted_pheno_diff", "drug_statin_TC_l_m_GWASformatted_pheno_log_diff",
            "TC_EUR_2021", "TC_EUR_2021",
            "drug_metformin_HbA1c_l_m_GWASformatted_pheno_diff", "drug_metformin_HbA1c_l_m_GWASformatted_pheno_log_diff",
            "drug_metformin_HbA1c_l_m_GWASformatted_pheno_diff", "drug_metformin_HbA1c_l_m_GWASformatted_pheno_log_diff",
            "HbA1c", "HbA1c",
            "drug_antihpt_all_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_all_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_all_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_all_SBP_l_m_GWASformatted_pheno_log_diff",
            "SBP_meta", "SBP_meta",
            "drug_antihpt_ACEi_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_ACEi_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_ACEi_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_ACEi_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_thiazide_diuretics_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_thiazide_diuretics_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_thiazide_diuretics_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_thiazide_diuretics_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_CCB_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_CCB_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_CCB_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_CCB_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_betablocker_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_betablocker_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_antihpt_betablocker_SBP_l_m_GWASformatted_pheno_diff", "drug_antihpt_betablocker_SBP_l_m_GWASformatted_pheno_log_diff",
            "drug_betablocker_HR_l_m_GWASformatted_pheno_diff", "drug_betablocker_HR_l_m_GWASformatted_pheno_log_diff",
            "drug_betablocker_HR_l_m_GWASformatted_pheno_diff", "drug_betablocker_HR_l_m_GWASformatted_pheno_log_diff",
            "PR", "PR"
            )
trait_names_1 = c("LDL-C response to statin (post-base)", "LDL-C response to statin (log(post)-log(base))",
                  "LDL-C response to statin (post-base)", "LDL-C response to statin (log(post)-log(base))",
                  "LDL-C (baseline)", "LDL-C (baseline)",
                  "HDL-C response to statin (post-base)", "HDL-C response to statin (log(post)-log(base))",
                  "HDL-C response to statin (post-base)", "HDL-C response to statin (log(post)-log(base))",
                  "HDL-C (baseline)", "HDL-C (baseline)",
                  "TC response to statin (post-base)", "TC response to statin (log(post)-log(base))",
                  "TC response to statin (post-base)", "TC response to statin (log(post)-log(base))",
                  "TC (baseline)", "TC (baseline)",
                  "HbA1c response to metformin (post-base)", "HbA1c response to metformin (log(post)-log(base))",
                  "HbA1c response to metformin (post-base)", "HbA1c response to metformin (log(post)-log(base))",
                  "HbA1c (baseline)", "HbA1c (baseline)",
                  "SBP response to antihpt (all) (post-base)", "SBP response to antihpt (all) (log(post)-log(base))",
                  "SBP response to antihpt (all) (post-base)", "SBP response to antihpt (all) (log(post)-log(base))",
                  "SBP (baseline)", "SBP (baseline)",
                  "SBP response to antihpt (ACEi) (post-base)", "SBP response to antihpt (ACEi) (log(post)-log(base))",
                  "SBP response to antihpt (ACEi) (post-base)", "SBP response to antihpt (ACEi) (log(post)-log(base))",
                  "SBP response to antihpt (diuretics) (post-base)", "SBP response to antihpt (diuretics) (log(post)-log(base))",
                  "SBP response to antihpt (diuretics) (post-base)", "SBP response to antihpt (diuretics) (log(post)-log(base))",
                  "SBP response to antihpt (CCB) (post-base)", "SBP response to antihpt (CCB) (log(post)-log(base))",
                  "SBP response to antihpt (CCB) (post-base)", "SBP response to antihpt (CCB) (log(post)-log(base))",
                  "SBP response to betablocker (post-base)", "SBP response to betablocker (log(post)-log(base))",
                  "SBP response to betablocker (post-base)", "SBP response to betablocker (log(post)-log(base))",
                  "HR response to betablocker (post-base)", "HR response to betablocker (log(post)-log(base))",
                  "HR response to betablocker (post-base)", "HR response to betablocker (log(post)-log(base))",
                  "HR (baseline)", "HR (baseline)"
)

traits_2 = c("LDL_EUR_2021", "LDL_EUR_2021",
            "control_statin_LDL_l_1_GWASformatted_pheno_diff", "control_statin_LDL_l_1_GWASformatted_pheno_log_diff",
            "control_statin_LDL_l_1_GWASformatted_pheno_diff", "control_statin_LDL_l_1_GWASformatted_pheno_log_diff",
            "HDL_EUR_2021", "HDL_EUR_2021",
            "control_statin_HDL_l_1_GWASformatted_pheno_diff", "control_statin_HDL_l_1_GWASformatted_pheno_log_diff",
            "control_statin_HDL_l_1_GWASformatted_pheno_diff", "control_statin_HDL_l_1_GWASformatted_pheno_log_diff",
            "TC_EUR_2021", "TC_EUR_2021",
            "control_statin_TC_l_1_GWASformatted_pheno_diff", "control_statin_TC_l_1_GWASformatted_pheno_log_diff",
            "control_statin_TC_l_1_GWASformatted_pheno_diff", "control_statin_TC_l_1_GWASformatted_pheno_log_diff",
            "HbA1c", "HbA1c",
            "control_metformin_HbA1c_l_1_GWASformatted_pheno_diff", "control_metformin_HbA1c_l_1_GWASformatted_pheno_log_diff",
            "control_metformin_HbA1c_l_1_GWASformatted_pheno_diff", "control_metformin_HbA1c_l_1_GWASformatted_pheno_log_diff",
            "SBP_meta", "SBP_meta",
            "control_antihpt_all_SBP_l_1_GWASformatted_pheno_diff", "control_antihpt_all_SBP_l_1_GWASformatted_pheno_log_diff",
            "control_antihpt_all_SBP_l_1_GWASformatted_pheno_diff", "control_antihpt_all_SBP_l_1_GWASformatted_pheno_log_diff",
            "SBP_meta", "SBP_meta",
            "control_antihpt_all_SBP_l_1_GWASformatted_pheno_diff", "control_antihpt_all_SBP_l_1_GWASformatted_pheno_log_diff",
            "SBP_meta", "SBP_meta",
            "control_antihpt_all_SBP_l_1_GWASformatted_pheno_diff", "control_antihpt_all_SBP_l_1_GWASformatted_pheno_log_diff",
            "SBP_meta", "SBP_meta",
            "control_antihpt_all_SBP_l_1_GWASformatted_pheno_diff", "control_antihpt_all_SBP_l_1_GWASformatted_pheno_log_diff",
            "SBP_meta", "SBP_meta",
            "control_antihpt_all_SBP_l_1_GWASformatted_pheno_diff", "control_antihpt_all_SBP_l_1_GWASformatted_pheno_log_diff",
            "PR", "PR",
            "control_betablocker_HR_l_1_GWASformatted_pheno_diff", "control_betablocker_HR_l_1_GWASformatted_pheno_log_diff",
            "control_betablocker_HR_l_1_GWASformatted_pheno_diff", "control_betablocker_HR_l_1_GWASformatted_pheno_log_diff"
            )

trait_names_2 = c("LDL-C (baseline)", "LDL-C (baseline)",
                  "LDL-C progression (post-base)", "LDL-C progression (log(post)-log(base))",
                  "LDL-C progression (post-base)", "LDL-C progression (log(post)-log(base))",
                  "HDL-C (baseline)", "HDL-C (baseline)",
                  "HDL-C progression (post-base)", "HDL-C progression (log(post)-log(base))",
                  "HDL-C progression (post-base)", "HDL-C progression (log(post)-log(base))",
                  "TC (baseline)", "TC (baseline)",
                  "TC progression (post-base)", "TC progression (log(post)-log(base))",
                  "TC progression (post-base)", "TC progression (log(post)-log(base))",
                  "HbA1c (baseline)", "HbA1c (baseline)",
                  "HbA1c progression (post-base)", "HbA1c progression (log(post)-log(base))",
                  "HbA1c progression (post-base)", "HbA1c progression (log(post)-log(base))",
                  "SBP (baseline)", "SBP (baseline)",
                  "SBP progression (post-base)", "SBP progression (log(post)-log(base))",
                  "SBP progression (post-base)", "SBP progression (log(post)-log(base))",
                  "SBP (baseline)", "SBP (baseline)",
                  "SBP progression (post-base)", "SBP progression (log(post)-log(base))",
                  "SBP (baseline)", "SBP (baseline)",
                  "SBP progression (post-base)", "SBP progression (log(post)-log(base))",
                  "SBP (baseline)", "SBP (baseline)",
                  "SBP progression (post-base)", "SBP progression (log(post)-log(base))",
                  "SBP (baseline)", "SBP (baseline)",
                  "SBP progression (post-base)", "SBP progression (log(post)-log(base))",
                  "HR (baseline)", "HR (baseline)",
                  "HR progression (post-base)", "HR progression (log(post)-log(base))",
                  "HR progression (post-base)", "HR progression (log(post)-log(base))"
)

colnames <- c("Trait1", "h2_1", "se_h2_1", "Trait2", "h2_2", "se_h2_2", "rg", "se_rg")
result_df <- data.frame(matrix(NA, length(traits_1), length(colnames)))
names(result_df) <- colnames

for (i in 1:length(traits_1)){
    print(paste0("Trait1: ", traits_1[i], "; Trait2: ", traits_2[i]))
    trait.names = c(traits_1[i], traits_2[i])
    traits = paste0(trait.names, ".sumstats.gz")
    sample.prev = rep(NA,length(traits)) #continuous traits
    population.prev = rep(NA,length(traits)) #continuous traits


    ldsc_result <- function(traits, sample.prev, population.prev, ld, trait.names){
        tryCatch(
            {
                invisible(utils::capture.output(LDSCoutput <- GenomicSEM::ldsc(traits,
                                                               sample.prev,
                                                               population.prev,
                                                               ld,
                                                               trait.names,
                                                               stand=TRUE)))
                return("ok")
            },
                warning = function(w) {
                return("w")
            }
        )
    }
    # run once, just to see if there is a warning
    warn = ldsc_result(traits, sample.prev, population.prev, ld, trait.names)

    # run to get the results
    invisible(utils::capture.output(LDSCoutput <- GenomicSEM::ldsc(traits,
                                                               sample.prev,
                                                               population.prev,
                                                               ld,
                                                               trait.names,
                                                               stand=TRUE)))

    S_Stand = LDSCoutput$S_Stand
    V_Stand = LDSCoutput$V_Stand

    r<-nrow(S_Stand)
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))

    # get heritabilities
    V = LDSCoutput$V
    S = LDSCoutput$S

    # get se
    SE<-matrix(0, r, r)
    SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(V))

    if (warn == "ok"){
        rg = S_Stand[2,1]
        se_rg = SE_Stand[2,1] 
    } else {
        rg = NA
        se_rg = NA
    }
                                                                  
    result_df[i, "Trait1"] = trait_names_1[i]
    result_df[i, "Trait2"] = trait_names_2[i]
    result_df[i, "h2_1"] = S[1,1]
    result_df[i, "h2_2"] = S[2,2]
    result_df[i, "se_h2_1"] = SE[1,1]
    result_df[i, "se_h2_2"] = SE[2,2]
    result_df[i, "rg"] = rg
    result_df[i, "se_rg"] = se_rg
}

#write.table(result_df, snakemake@output[[1]], sep = '\t', row.names = F, quote = F)
write.table(result_df, "genetic_correlations.tsv", sep = '\t', row.names = F, quote = F)

#### end ####
