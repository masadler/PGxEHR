import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input["file"], sep = '\t')
ref_snp_df = pd.read_csv(snakemake.input["ref_snps"], names = ["CHROM", "SNP", "d", "GENPOS", "A1", "A2"], sep = '\t')
ref_snps = list(ref_snp_df.SNP)

df = df.rename(columns = {"ID": "SNP", "ALLELE1": "A1", "ALLELE0": "A2", "A1FREQ": "Freq"})

# MAF filter
df = df[(df.Freq >= 0.05) & (df.Freq <= 0.95)]

# REF SNP filter
df = df[df.SNP.isin(ref_snps)]

df['p'] = 10**(-df.LOG10P)

print("Number of rows:", df.shape[0])

df = df[['SNP', 'CHROM', 'GENPOS', 'A1', 'A2', 'Freq', 'BETA', 'SE', 'p', 'N']]
df = df.dropna()

df[['SNP', 'CHROM', 'GENPOS', 'A1', 'A2', 'Freq', 'BETA', 'SE', 'p', 'N']].to_csv(snakemake.output[0], sep = '\t', index = False)
