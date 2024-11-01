library(dplyr)
library(stringr)

m = snakemake@wildcards[["measure"]]
filter = snakemake@wildcards[["filter"]]
n_meas = snakemake@wildcards[["n_meas"]]

df = read.csv(snakemake@input[["pheno"]], sep = '\t')

# define minimum baseline value per measure

if (m == "LDL"){
    l = 2
} else if (m == "HDL"){
    l = 0.2
} else if (m == "TC"){
    l = 5
}

# define stringent and lenient filtering specific covariates
if (filter == "s"){
    df$dose = df$dose_start
    covs_filter = ""
} else if (filter == "l"){
    df$dose = df$dose_avg
    covs_filter = "antilipemic_notstatin_prior + "
}

# filter on minimum baseline value
df = df[df$baseline_measure >= l,]

# determine which drugs to use, minimum 20 per type
med_count = df %>% count(drug_start)
drugs = as.vector(med_count[med_count$n >= 20, "drug_start"][[1]])
df = df[df$drug_start %in% drugs, ]

# define GWAS outcome phenotypes
df$pheno_diff_out = scale(df$post_measure - df$baseline_measure)
df$pheno_log_diff_out = scale(log(df$post_measure) - log(df$baseline_measure))
df$pheno_prop_out = scale(df$post_measure/df$baseline_measure)

# calculate residuals with relevant covariates

for (pheno in c("diff", "log_diff", "prop")){
    print(pheno)
    df[paste0("pheno_", pheno)] = residuals(lm(paste0("pheno_", pheno, "_out ~ sex + age_drug + drug_start*dose + post_measure_DT + ", 
                                    covs_filter, paste0("PC", 1:20, collapse="+")),
                                    data = df, na.action = na.exclude))

}

# classical literature PGx residual adjusted for baseline
pheno = "diff"
df[paste0("pheno_", pheno, "_adj")] = residuals(lm(paste0("pheno_", pheno, "_out ~ baseline_measure + sex + age_drug + drug_start*dose + post_measure_DT + ", 
                                    covs_filter, paste0("PC", 1:20, collapse="+")),
                                    data = df, na.action = na.exclude))

# format data for regenie and save
df$FID = df$eid
df$IID = df$eid

write.table(df[, c("FID", "IID", "pheno_diff", "pheno_log_diff", "pheno_prop", "pheno_diff_adj")], snakemake@output[[1]], row.names = F, quote = F, sep = " ", na = "NA")
