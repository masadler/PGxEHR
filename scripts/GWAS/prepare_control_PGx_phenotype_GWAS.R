library(dplyr)
library(stringr)

df = read.csv(snakemake@input[["pheno"]], sep = '\t')

# define GWAS outcome phenotypes
df$pheno_diff_out = scale(df$measure2 - df$measure1)
df$pheno_log_diff_out = scale(log(df$measure2) - log(df$measure1))
df$pheno_prop_out = scale(df$measure2/df$measure1)

# calculate residuals with relevant covariates

for (pheno in c("diff", "log_diff", "prop")){
    df[paste0("pheno_", pheno)] = residuals(lm(paste0("pheno_", pheno, "_out ~ sex + age_m1 + measure_DT + ", paste0("PC", 1:20, collapse="+")),
                                    data = df, na.action = na.exclude))
}

# classical literature PGx residual adjusted for baseline
pheno = "diff"
df[paste0("pheno_", pheno, "_adj")] = residuals(lm(paste0("pheno_", pheno, "_out ~ measure1 + sex + age_m1 + measure_DT + ", paste0("PC", 1:20, collapse="+")),
                                    data = df, na.action = na.exclude))

# format data for regenie and save
df$FID = df$eid
df$IID = df$eid

write.table(df[, c("FID", "IID", "pheno_diff", "pheno_log_diff", "pheno_prop", "pheno_diff_adj")], snakemake@output[[1]], row.names = F, quote = F, sep = " ", na = "NA")
