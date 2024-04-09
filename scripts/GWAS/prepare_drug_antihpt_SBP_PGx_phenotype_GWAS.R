library(dplyr)
library(stringr)

m = "SBP"
filter = snakemake@wildcards[["filter"]]
n_meas = snakemake@wildcards[["n_meas"]]
antihpt_class = snakemake@wildcards["medication"]

df = read.csv(snakemake@input[["pheno"]], sep = '\t')

# define minimum baseline value per measure

if (m == "SBP"){
    l = 120
}

# define stringent and lenient filtering specific covariates
if (filter == "s"){
    df$dose = df$dose_start
    if (antihpt_class == "all"){
        covs_filter = "drug_start + "
    } else {
        covs_filter = "drug_start*dose + "
    }
} else if (filter == "l"){
    df$dose = df$dose_avg
    if (antihpt_class == "all"){
        covs_filter = "drug_start + bblocker_prior + loop_diuretics_prior + "
    } else if (antihpt_class == "betablocker"){
        covs_filter = "drug_start*dose + loop_diuretics_prior + "
    } else {
        covs_filter = "drug_start*dose + bblocker_prior + loop_diuretics_prior + "
    }
}

# filter on minimum baseline value
df = df[df$baseline_measure >= l,]

# determine which drugs to use, minimum 20 per type
med_count = df %>% count(drug_start)
drugs = as.vector(med_count[med_count$n >= 20, "drug_start"])
df = df[df$drug_start %in% drugs, ]

# make sure there are enough levels for loop_diuretics_prior
if ((nrow(df[df$loop_diuretics_prior == "yes", ]) < 1) & (filter == "l")){
    covs_filter = unlist(strsplit(covs_filter, split= 'loop_diuretics_prior + ', fixed=TRUE))[1]
}

# define GWAS outcome phenotypes
df$pheno_diff_out = scale(df$post_measure - df$baseline_measure)
df$pheno_log_diff_out = scale(log(df$post_measure) - log(df$baseline_measure))
df$pheno_prop_out = scale(df$post_measure/df$baseline_measure)

# calculate residuals with relevant covariates

for (pheno in c("diff", "log_diff", "prop")){
    df[paste0("pheno_", pheno)] = residuals(lm(paste0("pheno_", pheno, "_out ~ BMI*sex + age_drug + post_measure_DT + ", 
                                    covs_filter, paste0("PC", 1:20, collapse="+")),
                                    data = df, na.action = na.exclude))
}

# classical literature PGx residual adjusted for baseline
pheno = "diff"
df[paste0("pheno_", pheno, "_adj")] = residuals(lm(paste0("pheno_", pheno, "_out ~ baseline_measure + BMI*sex + age_drug + post_measure_DT + ", 
                                    covs_filter, paste0("PC", 1:20, collapse="+")),
                                    data = df, na.action = na.exclude))

# format data for regenie and save
df$FID = df$eid
df$IID = df$eid

write.table(df[, c("FID", "IID", "pheno_diff", "pheno_log_diff", "pheno_prop", "pheno_diff_adj")], snakemake@output[[1]], row.names = F, quote = F, sep = " ", na = "NA")
