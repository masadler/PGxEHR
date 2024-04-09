library(rbgen)
library(abind)

# read sample information
df = read.csv(snakemake@input[["cov"]], sep = '\t')

# map the people to order
sample_df = read.table(snakemake@input[["sample_order"]], header = T)
sample_df = sample_df[2:nrow(sample_df), ]# delete second header

## add column with sample name from bgen files
sample_df["order"] = as.character(seq.int(nrow(sample_df)))
sample_df["bgen_ID"] = paste0("(anonymous_sample_", sample_df$order, ")")

## samples of interest
samples = as.vector(sample_df[sample_df$ID_1 %in% df$eid, "bgen_ID"])
print(paste0("Number of samples: ", length(samples)))

rsid = snakemake@wildcards[["snp"]]
file = bgen.load(snakemake@input[["geno"]], rsids = rsid, samples = samples)

genotype_df = file$data[,,c(2,3), drop = F] # we don't need g0
variants_df = file$variants

genotype_array = aperm(as.array(genotype_df[rsid,,,drop = F]), c(2,1,3))
genotype_array = abind(genotype_array, genotype_array[,,1,drop = F]+2*genotype_array[,,2,drop = F], along = 3) # get gene dosage for each person, 0, 1, 2 encoding

dosage_mat = as.matrix(genotype_array[,,3])
dosage_mat[is.nan(dosage_mat)] = 0

result_df = data.frame(bgen_ID = samples, dosage = dosage_mat[, 1])
result_df = merge(result_df, sample_df[, c("ID_1", "bgen_ID")], by = "bgen_ID")
result_df$A1 = variants_df[1, "allele1"]
result_df$A2 = variants_df[1, "allele0"]

result_df = result_df[, c("ID_1", "A1", "A2", "dosage")]
names(result_df) = c("eid", "A1", "A2", "dosage")

write.table(result_df, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')





