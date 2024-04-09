import pandas as pd

result_df = pd.DataFrame()

for file in snakemake.input["files"]:
    df = pd.read_csv(file, delim_whitespace=True)
    result_df = pd.concat([result_df, df])

result_df.to_csv(snakemake.output[0], sep = '\t', index = False)