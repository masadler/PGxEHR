library(dplyr)
library(stringr)

heatmap_res <- function(df, x, y, z){

    df <- df %>%
    mutate(
        base_quintile = ntile(!!sym(x), 5),
        base_quintile = as.factor(base_quintile)
    )

    # first by baseline and then by PRS

    result_mat = matrix(NA, 5, 5)

    for (i in 1:5){
        df_quin = df[df$base_quintile == i, ]
        df_quin <- df_quin %>%
                    mutate(
                        PRS_quintile = ntile(!!sym(y), 5),
                        PRS_quintile = as.factor(PRS_quintile)
                    )
        sum_df <- df_quin %>% 
                group_by(PRS_quintile) %>% 
                summarise(mean = mean(!!sym(z)),
                          median = median(!!sym(z)),
                          sd = sd(!!sym(z)))

        result_mat[, i] = as.vector(data.frame(sum_df[, "median"]))[[1]]
        #result_mat[, i] = as.vector(data.frame(sum_df[, "mean"]))[[1]]
    }

    result_df = data.frame(result_mat)
    names(result_df) = paste0(seq(1,5))
    row.names(result_df) = paste0(seq(1,5))
    return(result_df)
}

#### statin - LDL: LDL PRS ####

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

# filter on minimum baseline value
df = df[df$baseline_measure >= l,]

# determine which drugs to use, minimum 20 per type
med_count = df %>% count(drug_start)
drugs = as.vector(med_count[med_count$n >= 20, "drug_start"])
df = df[df$drug_start %in% drugs, ]

# add PRS of measure
prs_df = read.csv(snakemake@input[["prs"]], sep = '\t')
prs_trait_df = prs_df[, c("IID", paste0(m, "_UKBB"))]
names(prs_trait_df) = c("eid", "PRS")
prs_trait_df = prs_trait_df[!is.na(prs_trait_df$PRS), ]

df = merge(df, prs_trait_df, by = "eid")

# define PGx drug response
df$pheno_diff = df$post_measure - df$baseline_measure
df$pheno_prop = df$post_measure/df$baseline_measure
df$pheno_prop_log = log(df$post_measure/df$baseline_measure)

# calculate baseline adjusted for PRS
df$baseline_prs_res = residuals(lm("baseline_measure ~ PRS", data = df, na.action = na.exclude))

# calculate heatmap for diff drug response
result_df = heatmap_res(df, "baseline_prs_res", "PRS", "pheno_diff")
write.table(result_df, snakemake@output[["absdiff_baseres_heatmap"]], row.names = F, quote = F, sep = '\t')

result_df = heatmap_res(df, "baseline_measure", "PRS", "pheno_diff")
write.table(result_df, snakemake@output[["absdiff_base_heatmap"]], row.names = F, quote = F, sep = '\t')

# calculate heatmap for ratio (proportion) drug response
result_df = heatmap_res(df, "baseline_prs_res", "PRS", "pheno_prop")
write.table(result_df, snakemake@output[["ratio_baseres_heatmap"]], row.names = F, quote = F, sep = '\t')

result_df = heatmap_res(df, "baseline_measure", "PRS", "pheno_prop")
write.table(result_df, snakemake@output[["ratio_base_heatmap"]], row.names = F, quote = F, sep = '\t')

# calculate heatmap for logarithmic difference drug response
result_df = heatmap_res(df, "baseline_prs_res", "PRS", "pheno_prop_log")
write.table(result_df, snakemake@output[["logdiff_baseres_heatmap"]], row.names = F, quote = F, sep = '\t')

result_df = heatmap_res(df, "baseline_measure", "PRS", "pheno_prop_log")
write.table(result_df, snakemake@output[["logdiff_base_heatmap"]], row.names = F, quote = F, sep = '\t')


