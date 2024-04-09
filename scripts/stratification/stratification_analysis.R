library(dplyr)
library(data.table)
library(stringr)

colnames = c("trait", "adjusted", "SNP", "gene", "time", "drug", "genotype", "mean", "sd", "N")
result_df = data.frame(matrix(NA, nrow = 0, ncol = length(colnames)))
names(result_df) = colnames

colnames = c("trait", "adjusted", "SNP", "gene", "group", "b", "se", "p", "N")
reg_df = data.frame(matrix(NA, nrow = 0, ncol = length(colnames)))
names(reg_df) = colnames

# stratification analysis
medication_measure_pairs = c("statin_LDL", "statin_LDL", "statin_TC", "statin_HDL")
snps = c("rs7412", "rs4970837", "rs4149056", "rs11076175")
genes = c("APOE", "SORT1", "SLCO1B1", "CETP")
allele1 = c("T", "T", "C", "G")
allele2 = c("C", "G", "T", "A")

i = 1
r= 1

for (pair in medication_measure_pairs){
    print(pair)

    # prepare drug data
    df = read.csv(paste0("output/PGx_phenotype/drug_", pair, "_l_m_PGx_phenotype.tsv"), sep = '\t')

    if (pair == "statin_LDL"){
        unit = "mmol/L"
        l = 2
        pair_name = "LDL-C response to statin"
    } else if (pair == "statin_TC"){
        unit = "mmol/L"
        l = 5
        pair_name = "TC response to statin"
    } else if (pair == "statin_HDL"){
        unit = "mmol/L"
        l = 0.2
        pair_name = "HDL-C response to statin"
    } 

    df = df[df$baseline_measure > l,]

    #med_count = df %>% count(drug_start)
    #drugs = as.vector(med_count[med_count$n >= 20, "drug_start"])
    #df = df[df$drug_start %in% drugs, ]
    df = df[(df$drug_start == "simvastatin") & (df$dose_start == 40), ]
    df$drug = 1

    # prepare control data
    control_df = read.csv(paste0("output/PGx_phenotype/control_", pair, "_l_1_PGx_phenotype.tsv"), sep = '\t')
    control_df = control_df %>% rename(post_measure = measure2, baseline_measure = measure1, age_drug = age_m1)
    control_df$drug = 0
    
    # combine the two
    cols = c("eid", "post_measure", "baseline_measure", "sex", "age_drug", "drug")

    data_df = rbind(df[, cols], control_df[, cols])
    data_df$pheno_diff = data_df$post_measure - data_df$baseline_measure
    data_df$pheno_log_diff = log(data_df$post_measure) - log(data_df$baseline_measure)

    for (a in c(0,1)){ # a = 0: unadjusted for baseline, a = 1: adjusted for baseline

        snp_df = fread(paste0("output/stratification/", snps[i], "_dosage.tsv"))

        df.snp = merge(data_df, snp_df)

        # stratification
        alleles = c(paste0(rep(allele1[i], 2), collapse = ""), paste0(allele1[i], allele2[i], collapse = ""), paste0(rep(allele2[i], 2), collapse = ""))

        for (drug in c(0,1,2)){ # drug = 0: control; drug = 1: drug; drug = 2: baseline
            
            if (drug == 2){
                time = "base"

                df.snp.group = df.snp

                # regression for baseline
                df.snp.group$baseline_measure_std = scale(log(df.snp.group$baseline_measure))

                df.snp.group["pheno"] <- residuals(lm("baseline_measure_std ~ sex + age_drug",
                                                        data = df.snp.group, na.action = na.exclude))

                reg <- summary(lm(baseline_measure_std ~ sex + age_drug + dosage, data = df.snp.group, na.action = na.exclude))
                
                reg_df[r, "trait"] = pair
                reg_df[r, "adjusted"] = "NA"
                reg_df[r, "SNP"] = snps[i]
                reg_df[r, "gene"] = genes[i]
                reg_df[r, "group"] = "baseline"
                reg_df[r, "b"] = reg$coef["dosage", "Estimate"]
                reg_df[r, "se"] = reg$coef["dosage", "Std. Error"]
                reg_df[r, "p"] = reg$coef["dosage", "Pr(>|t|)"]
                reg_df[r, "N"] = nrow(df.snp.group)
            } else {

                time = "post"

                # regression in each group
                df.snp.group = df.snp[df.snp$drug == drug, ]
                df.snp.group$log_diff_std = scale(df.snp.group$pheno_log_diff)

                if (a == 0){
                    df.snp.group["pheno"] <- residuals(lm("log_diff_std ~ sex + age_drug",
                                                            data = df.snp.group, na.action = na.exclude))
                    reg <- summary(lm(log_diff_std ~ sex + age_drug + dosage, data = df.snp.group, na.action = na.exclude))
                } else {
                    df.snp.group["pheno"] <- residuals(lm("log_diff_std ~ baseline_measure + sex + age_drug",
                                        data = df.snp.group, na.action = na.exclude))
                    reg <- summary(lm(log_diff_std ~ baseline_measure + sex + age_drug + dosage, data = df.snp.group, na.action = na.exclude))
                }

                reg_df[r, "trait"] = pair
                reg_df[r, "adjusted"] = a
                reg_df[r, "SNP"] = snps[i]
                reg_df[r, "gene"] = genes[i]
                reg_df[r, "group"] = drug
                reg_df[r, "b"] = reg$coef["dosage", "Estimate"]
                reg_df[r, "se"] = reg$coef["dosage", "Std. Error"]
                reg_df[r, "p"] = reg$coef["dosage", "Pr(>|t|)"]
                reg_df[r, "N"] = nrow(df.snp.group)
            }
            r = r+1

            pheno = "pheno"
            risk_category = list()
            risk_category[[1]] = df.snp.group[(df.snp.group$dosage >= 1.5), pheno]  # a1, a1
            risk_category[[2]] = df.snp.group[(df.snp.group$dosage < 1.5) & (df.snp.group$dosage > 0.5), pheno]  # a1, a2
            risk_category[[3]] = df.snp.group[(df.snp.group$dosage <= 0.5), pheno]  # a2, a2

            for (k in 1:3){
                row = i*18-18 + a*9 + drug*3 + k
                result_df[row, "trait"] = pair
                result_df[row, "adjusted"] = a
                result_df[row, "SNP"] = snps[i]
                result_df[row, "gene"] = genes[i]
                result_df[row, "time"] = time
                result_df[row, "drug"] = drug
                result_df[row, "genotype"] = alleles[k]
                result_df[row, "mean"] = mean(risk_category[[k]])
                result_df[row, "sd"] = sd(risk_category[[k]])
                result_df[row, "N"] = length(risk_category[[k]])
            }
        }
    }
    i = i+1
}

write.table(result_df, snakemake@output[[1]], sep = '\t', quote = F, row.names = F)
write.table(reg_df, snakemake@output[[2]], sep = '\t', quote = F, row.names = F)

