library(dplyr)
library(stringr)
library(data.table)

df = read.csv(snakemake@input[["pheno"]], sep = '\t')
prs_df = read.csv(snakemake@input[["prs_ukbb"]], sep = '\t')
snp_df = fread(snakemake@input[["snp"]])

prs_trait_df = prs_df[, c("IID", "LDL_UKBB")]
names(prs_trait_df) = c("eid", "PRS")
prs_trait_df = prs_trait_df[!is.na(prs_trait_df$PRS), ]

df = merge(df, prs_trait_df, by = "eid")

# filtering
l = 2
df = df[df$baseline_measure > l,]
med_count = df %>% count(drug_start)
drugs = as.vector(med_count[med_count$n >= 20, "drug_start"])
df = df[df$drug_start %in% drugs, ]

# add snp
df = merge(df, snp_df)

colnames = c("Baseline", "PRS", "dosage", "post_measure")
result_df = data.frame(matrix(NA, nrow = 0, ncol = length(colnames)))

df <- df %>%
    mutate(
        base_tertile = ntile(`baseline_measure`, 3)
    )

# first by baseline and then by PRS and then by APOE
for (i in 1:3){
    df_tert = df[df$base_tertile == i, ]
    df_tert <- df_tert %>%
                mutate(
                    PRS = ntile(`PRS`, 3),
                )

    df_tert$Baseline = i
    result_df = rbind(result_df, df_tert[, colnames])
}

write.table(result_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')

# numerical statistics

df = result_df

df = df[df$dosage != 2, ]
df$APOE = "CC"
df[df$dosage == 1, "APOE"] = "TC"

df$Baseline_name = "Baseline tertile 1"
df[df$Baseline == 2, "Baseline_name"] = "Baseline tertile 2"
df[df$Baseline == 3, "Baseline_name"] = "Baseline tertile 3"

df$Baseline_name = factor(df$Baseline_name)
df$PRS = factor(df$PRS)
df$APOE = factor(df$APOE, levels = c("TC", "CC"))

df$PRS_name = "PRS tertile 1"
df[df$PRS == 2, "PRS_name"] = "PRS tertile 2"
df[df$PRS == 3, "PRS_name"] = "PRS tertile 3"

summary_stats <- df %>%
  group_by(Baseline_name, PRS_name, APOE) %>%
  summarise(
    mean = mean(post_measure),
    sd = sd(post_measure),
    median = median(post_measure),
    N = n()
  )

write.table(summary_stats, snakemake@output[[2]], sep = '\t', quote = F, row.names = F)
