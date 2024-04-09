import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timedelta, date

"""
Control two timepoint HbA1c GWAS using individuals with no antidiabetic entries

2 measures: 6 months - 3 years apart

"""

m = snakemake.wildcards["measure"] # HbA1c
filter = snakemake.wildcards["filter"] # lenient only (in terms of time interval)
n_meas = snakemake.wildcards["n_meas"] # 1, two time points only

#### biomarker measures filtering parameters

codes = ["42W4.", "42W5.", "XaERp", "XaPbt", "UKBB"] # "42W4." and "XaERp" are already converted to mmol/mol
min = 20
max = 195

#### time between first and second measure ####
interval_time_min = 180
interval_time_max = 1095

#### filter biomarker measures
measure_df = pd.read_csv(snakemake.input["measure"], sep = '\t')
measure_df["event_dt"] = pd.to_datetime(measure_df["event_dt"], format = "%Y-%m-%d").dt.date
measure_df = measure_df[measure_df.read_code.isin(codes)]
measure_df = measure_df[(measure_df[m] >= min) & (measure_df[m] <= max)]
measure_df = measure_df.dropna(subset=[m, "event_dt"])
# drop duplicated measures on the same day
measure_df = measure_df.drop_duplicates(subset = ["eid", "event_dt"])

# eliminate individuals with antidiabetic entries
medication_df = pd.read_csv(snakemake.input["medication"], sep = '\t')
measure_df = measure_df[~measure_df.eid.isin(medication_df.eid)]

#### extract PGx control phenotype ####

keys = ["eid", "measure1_date", "measure2_date", "measure1", "measure2", "measure_DT"]
result = {key: [] for key in keys}

# document why individuals are removed
removed = {key: [] for key in ["eid", "reason"]}

for eid in measure_df.eid.unique():

    measure_person_df = measure_df[measure_df.eid == eid]
    
    if measure_person_df.shape[0] < 2:
        removed["eid"].append(eid)
        removed["reason"].append("no 2 measures")
        continue

    measure_person_df = measure_person_df.sort_values("event_dt")
    measure_person_df = measure_person_df.reset_index(drop = True)
    first_measure = measure_person_df.event_dt.iloc[0]
    measure_person_df["time"] = measure_person_df["event_dt"].apply(lambda x: (x - first_measure).days)

    two_m_interval_df = measure_person_df[(measure_person_df.time >= interval_time_min) &
                                          (measure_person_df.time <= interval_time_max)]
    
    if two_m_interval_df.shape[0] < 1:
        removed["eid"].append(eid)
        removed["reason"].append("no matching interval")
        continue
    
    # take the first set of measures fulfilling condition
    second_index = two_m_interval_df.index[0]
        
    # store result  
    result["eid"].append(eid)
    result["measure1_date"].append(measure_person_df["event_dt"].iloc[second_index-1])
    result["measure2_date"].append(measure_person_df["event_dt"].iloc[second_index])
    result["measure1"].append(measure_person_df[m].iloc[second_index-1])
    result["measure2"].append(measure_person_df[m].iloc[second_index])
    result["measure_DT"].append((measure_person_df["event_dt"].iloc[second_index]-measure_person_df["event_dt"].iloc[second_index-1]).days)

result_df = pd.DataFrame(result)
removed_df = pd.DataFrame(removed)

# add covariates
cov_df = pd.read_csv(snakemake.input["cov"], sep = '\t')
result_df = result_df.merge(cov_df, on = "eid")
result_df["day"] = 1 # dummy birth day
result_df = result_df.rename(columns = {"birth_year": "year", "birth_month": "month"})
result_df["birth_date"] = pd.to_datetime(result_df[["year", "month", "day"]])
result_df["age_m1"] = ((pd.to_datetime(result_df["measure1_date"]) - pd.to_datetime(result_df["birth_date"])).dt.days/365.25).astype(int)

result_df.to_csv(snakemake.output[0], sep = '\t', index = False)
removed_df.to_csv(snakemake.output[1], sep = '\t', index = False)