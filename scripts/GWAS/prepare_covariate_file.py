import numpy as np
import pandas as pd

cov_df = pd.read_csv(snakemake.input["cov"])

# remove withdrawn individuals
rm_df = pd.read_csv(snakemake.input["withdrawn"], names = ["eid"])
cov_df = cov_df[~cov_df.eid.isin(rm_df.eid)]

# included white unrelated
in_df = pd.read_csv(snakemake.input["included"], names = ["eid"])
cov_df = cov_df[cov_df.eid.isin(in_df.eid)]

#### extract phenotypes ####

# BMI
bmi_id = snakemake.params["bmi"]
year_id = snakemake.params["year"]
month_id = snakemake.params["month"]
usecols = ['eid'] + [bmi_id] + [year_id] + [month_id]
ukbb_df = pd.read_csv(snakemake.input["file"], usecols = usecols)
ukbb_df = ukbb_df.rename(columns = {bmi_id: "BMI", year_id: "year", month_id: "month"})

cov_df = cov_df.merge(ukbb_df, on = "eid")

cov_df.to_csv(snakemake.output[0], index = False, sep = '\t')

