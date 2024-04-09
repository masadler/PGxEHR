library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

sample_qc <- as.data.frame(fread(snakemake@input[["ukb_samples"]], header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

sample_eid <- as.data.frame(fread(snakemake@input[["sample_order"]], header = F, select = c(1,5), col.names = c("eid", "sex")))

df <- cbind(sample_eid, sample_qc)

#################################################
### Sample Filtering ############################

#################################################
### STEP 1: Exclude related samples (pca_calculation = 1)
#df <- df[which(df$pca_calculation == 1), ]
#print(paste0("STEP 1: Exclude related samples: ", nrow(df), " individuals"))

#################################################
### STEP 2: Exclude non-white, non-British ancestry samples (white_british = 1)
df <- df[which(df$white_british == 1), ]  
print(paste0("STEP 2: Exclude non-white, non-British ancestry samples: ", nrow(df), " individuals"))

### exclude samples with excess relatives
df <- df[which(df$excess_relatives == 0), ]  

#################################################
### STEP 3: Exclude redacted samples
df <- df[-which(df$eid < 0), ]  

#################################################
### STEP 7: Exclude samples with non-matching submitted vs. inferred sex
df <- df[which(df$submitted_gender == df$inferred_gender), ]
print(paste0("STEP 7: Exclude sex mismatches: ", nrow(df), " individuals"))

### STEP 8: Exclude withdrawn samples
withdrawn_df = read.table(snakemake@input[["withdrawn"]])
names(withdrawn_df) = c("eid")
df = df[!(df$eid %in% withdrawn_df$eid),]

write.table(df[, c("eid")], snakemake@output[[1]], col.names = F, quote = F, row.names = F)
