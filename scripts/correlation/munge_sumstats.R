library(GenomicSEM)
library(stringr)

#### munge summary statistics ####

snps = snakemake@input[["snps"]]

# baseline GWAS (GLGC)
glgc_gwas = snakemake@input[["glgc_gwas"]]

for (gwas in glgc_gwas){
    print(gwas)
    trait = str_split(tail(str_split(gwas, "/")[[1]], n = 1), "_summary_stats.tsv")[[1]][1]
    print(trait)
    #munge(files=c(gwas),hm3=snps, column.names=list(SNP = c("rs_id")), trait.names=c(paste0("output/ldsc/", trait)))
}

# baseline GWAS
base_gwas = snakemake@input[["base_gwas"]]

for (gwas in base_gwas){
    print(gwas)
    trait = str_split(tail(str_split(gwas, "/")[[1]], n = 1), "_gwas_summary_uk10kck.ma")[[1]][1]
    munge(files=c(gwas),hm3=snps,trait.names=c(paste0("output/ldsc/", trait)))
}

# PGx GWAS
pgx_gwas = snakemake@input[["pgx_gwas"]]

for (gwas in pgx_gwas){
    print(gwas)
    trait = str_split(tail(str_split(gwas, "/")[[1]], n = 1), ".tsv")[[1]][1]
    munge(files=c(gwas),hm3=snps,trait.names=c(paste0("output/ldsc/", trait)))
}

# longitudinal control GWAS
control_gwas = snakemake@input[["control_gwas"]]

for (gwas in control_gwas){
    print(gwas)
    trait = str_split(tail(str_split(gwas, "/")[[1]], n = 1), ".tsv")[[1]][1]
    munge(files=c(gwas),hm3=snps,trait.names=c(paste0("output/ldsc/", trait)))
}
#### end ####