import pandas as pd

clinical_data = pd.read_csv(snakemake.input["gp_clinical"], sep="\t", usecols = ["eid", "event_dt"], encoding= 'unicode_escape')
script_data = pd.read_csv(snakemake.input["gp_scripts"], sep="\t", usecols = ["eid", "issue_date"], encoding= 'unicode_escape')
reg_data = pd.read_csv(snakemake.input["gp_reg"], sep="\t", usecols = ["eid", "reg_date"], encoding= 'unicode_escape')

# get all primary care users
primary_data_df = pd.concat([clinical_data[["eid"]], script_data[["eid"]], reg_data[["eid"]]])
primary_data_df = primary_data_df.drop_duplicates()

primary_data_df.to_csv(snakemake.output[0], index = False, header = False)