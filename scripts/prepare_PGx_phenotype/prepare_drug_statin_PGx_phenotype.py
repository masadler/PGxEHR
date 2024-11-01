import numpy as np
import pandas as pd
from datetime import datetime
from datetime import timedelta, date

"""
Extract statin drug response phenotypes

"""

m = snakemake.wildcards["measure"]
filter = snakemake.wildcards["filter"]
n_meas = snakemake.wildcards["n_meas"]

#### biomarker measures filtering parameters

if m == "LDL":
    min = 0.26
    max = 10.3
    codes = ["44P6.", "44PI.", "XaIp4", "44PD.", "UKBB"]
elif m == "HDL":
    min = 0.2
    max = 3
    codes = ["44P5.", "44PB.", "UKBB"]
elif m == "TC":
    min = 1
    max = 15
    codes = ["44P..", "44PH.", "44PJ.", "XaJe9", "XE2eD", "UKBB"]

#### filtering scenario parameters; s-stringent, l-lenient

if filter == "s":

    # prescription regularity
    completeness = 0.6 # regular prescriptions 60% of time
    
    # time between baseline measures & prescription
    baseline_dt_before = 100
    baseline_dt_after = 7

    # time between prescription & post-measure
    post_dt_min = 100
    post_dt_max = 365

elif filter == "l":

    # prescription regularity
    completeness = 0.3 # regular prescriptions 30% of time
    
    # time between baseline measures & prescription
    baseline_dt_before = 365
    baseline_dt_after = 7

    # time between prescription & post-measure
    post_dt_min = 100
    post_dt_max = 730

#### other parameters
ti = 60 # interval in days between time points
prim_care_time = 730 # time before which an entry in primary care is necessary
dt_no_main = 365 # time before which no lipid-lowering prescription should be present

#### filter biomarker measures
measure_df = pd.read_csv(snakemake.input["measure"], sep = '\t')
measure_df["event_dt"] = pd.to_datetime(measure_df["event_dt"], format = "%Y-%m-%d").dt.date
measure_df = measure_df[measure_df.read_code.isin(codes)]
measure_df = measure_df[(measure_df[m] >= min) & (measure_df[m] <= max)]
measure_df = measure_df.dropna(subset=[m])
# drop duplicated measures on the same day
measure_df = measure_df.drop_duplicates(subset = ["eid", "event_dt"])

#### read medication data
medication_df = pd.read_csv(snakemake.input["medication"], sep = '\t')
medication_df["issue_date"] = pd.to_datetime(medication_df["issue_date"]).dt.date
medication_df = medication_df.replace('nan', np.NaN)
medication_df.combi = medication_df.combi.fillna('no')

# impute dose with median of statin type if missing (7 out of 4,165,780 times it is missing for statins)
#statins = ["atorvastatin", "simvastatin", "rosuvastatin", "pravastatin", "fluvastatin", "cerivastatin", "lovastatin", "pitavastatin"]
medication_df["is_statin"] = ~medication_df.drug.isna()
median_dose_by_type = medication_df[medication_df.is_statin == True].groupby('drug')['dose'].median()
medication_df = medication_df.rename(columns = {"dose": "dose_original"})
medication_df['dose'] = medication_df.apply(lambda row: 
                                               median_dose_by_type[row['drug']] 
                                               if (pd.isnull(row['dose_original']) & row.is_statin & (row['combi'] == "no"))
                                               else row['dose_original']
                                               , axis=1)

#### read all primary data for QC1 (1 entry 2 years before)
clinical_data = pd.read_csv(snakemake.input["gp_clinical"], sep="\t", usecols = ["eid", "event_dt"], encoding= 'unicode_escape')
clinical_data["date"] = pd.to_datetime(clinical_data["event_dt"], format = "%d/%m/%Y").dt.date

script_data = pd.read_csv(snakemake.input["gp_scripts"], sep="\t", usecols = ["eid", "issue_date"], encoding= 'unicode_escape')
script_data["date"] = pd.to_datetime(script_data["issue_date"], format = "%d/%m/%Y").dt.date

# combine data
primary_data_df = pd.concat([clinical_data[["eid", "date"]], script_data[["eid", "date"]]])
primary_data_df = primary_data_df.drop_duplicates()

#### extract PGx phenotype ####

def med_qc(df, time, drug, dose, post_times, filter):
    
    qc = True
    compl = 0
    reason = ""
    first_post_time = post_times[0]
    max_post_time = first_post_time
    dose_avg = dose

    # check all available post-measurement times, break when there is a change in medication prescriptions
    for post_time in post_times:
        # prescription available after post-measure (part of QC4)
        df_post = df[df.issue_date > post_time]
        if df_post.shape[0] < 1:
            if post_time != first_post_time:
                break
            qc = False
            reason = "no prescription after post-measure"
            return qc, compl, dose_avg, reason, max_post_time

        # prescriptions between initiation and post-measure
        df_during = df[df.issue_date < post_time]

        # QC5: no medication changes (dose excluded)
        if sum(df_during["combi"] != "no") > 0:
            if post_time != first_post_time:
                break
            qc = False
            reason = "combination therapy"
            return qc, compl, dose_avg, reason, max_post_time

        if filter == "s":
            df_same = df_during[(df_during.drug == drug) & (df_during.dose == dose)]
        else:
            df_same = df_during[(df_during.drug == drug)]
            
        if df_during.shape[0] != df_same.shape[0]:
            if post_time != first_post_time:
                break
            qc = False
            reason = "drug/dose change"
            return qc, compl, dose_avg, reason, max_post_time
        
        max_post_time = post_time
        if qc == False: # only relevant when post_time is the first post-measurement
            break

    df = df[df.issue_date < max_post_time]
    post_time = max_post_time

    # calculate average dose (can be the same than starting dose)
    if filter == "l":
        dose_avg = np.mean(df.dose)

    # QC4: regularity
    if qc == True:
        
        df["time"] = df["issue_date"].apply(lambda x: (x - time).days)
        tp = int((post_time-time).days / ti) # grid of every 2 months until post-measure
        prescription_grid = []
        for t in range(1, (tp+1)):
            df_int = df[(df.time >= (t*ti - ti/2)) & (df.time < (t*ti + ti/2))]
            if df_int.shape[0] == 0:
                prescription_grid.append(0)
            else:
                prescription_grid.append(1)
        compl = sum(prescription_grid)/len(prescription_grid)
        if compl < completeness:
            qc = False
            reason = "prescription missingness too high"

    return qc, compl, dose_avg, reason, max_post_time

def PGx_phenotype(medication_df, measure_df, primary_data_df, filter, n_meas, baseline_dt_before, baseline_dt_after, post_dt_min, post_dt_max, prim_care_time):

    keys = ["eid", "drug_start_date", "drug_start", "dose_start", "dose_avg", "baseline_measure", "baseline_DT", "post_measure", "post_measure_DT", "completeness", "antilipemic_notstatin_prior"]

    result = {key: [] for key in keys}

    # document why individuals are removed
    removed = {key: [] for key in ["eid", "reason"]}

    for eid in medication_df.eid.unique():

        antilipemic_notstatin_prior = "no"

        medication_person_df = medication_df[medication_df.eid == eid]
        measure_person_df = measure_df[measure_df.eid == eid]
        primary_data_person_df = primary_data_df[primary_data_df.eid == eid]
        
        if measure_person_df.shape[0] < 1:
            removed["eid"].append(eid)
            removed["reason"].append("no phenotype measure within range")
            continue
        
        # medication start - essential & secondary
        medication_person_df_essential = medication_person_df[~medication_person_df.drug.isna()]

        if medication_person_df_essential.shape[0] < 1: # no statin prescription at all
            removed["eid"].append(eid)
            removed["reason"].append("no primary medication")
            continue

        medication_person_df_nonessential = medication_person_df[medication_person_df.drug.isna()]

        medication_person_df_essential = medication_person_df_essential.sort_values("issue_date")
        start_time = medication_person_df_essential.issue_date.iloc[0]
        drug_start_info = medication_person_df_essential.iloc[0]

        # QC3
        if drug_start_info["combi"] != "no":
            removed["eid"].append(eid)
            removed["reason"].append("combination medication")
            continue

        # QC3
        if medication_person_df_nonessential.shape[0] > 0:
            medication_person_df_nonessential = medication_person_df_nonessential.sort_values("issue_date")

            if medication_person_df_nonessential[(medication_person_df_nonessential.issue_date <= (start_time)) & 
                                                (medication_person_df_nonessential.issue_date  > (start_time - timedelta(days=dt_no_main)))].shape[0] > 0:

                antilipemic_notstatin_prior = "yes"

                if filter == "s":
                    removed["eid"].append(eid)
                    removed["reason"].append("prior antilipemic use")
                    continue

        medication_person_df = medication_person_df[medication_person_df.issue_date > start_time] # successive medication
        medication_person_df = medication_person_df.sort_values("issue_date")

        # QC2 - baseline measures (preceeding medication start)
        baseline_measure = measure_person_df[(measure_person_df.event_dt <= (start_time + timedelta(days=baseline_dt_after))) &
                                            (measure_person_df.event_dt >= (start_time - timedelta(days=baseline_dt_before)))]
        baseline_measure = baseline_measure.sort_values("event_dt")
        
        if baseline_measure.shape[0] < 1:
            removed["eid"].append(eid)
            removed["reason"].append("baseline-medication time")
            continue
        
        baseline_time = baseline_measure.event_dt.iloc[-1] # take last measure as baseline time

        # QC1
        min_time = start_time
        if start_time < baseline_time:
            min_time = baseline_time
        
        # any entry should be present in primary care data up to a certain time before
        prim_care = primary_data_person_df[(primary_data_person_df.date <= (min_time - timedelta(days=prim_care_time)))]
        
        if prim_care.shape[0] < 1:
            removed["eid"].append(eid)
            removed["reason"].append("no prior primary care entry")
            continue
        
        # post measures within defined window
        post_measure = measure_person_df[(measure_person_df.event_dt >= (start_time + timedelta(days=post_dt_min))) &
                                        (measure_person_df.event_dt <= (start_time + timedelta(days=post_dt_max)))]
        post_measure = post_measure.sort_values("event_dt")
        
        if post_measure.shape[0] < 1:
            removed["eid"].append(eid)
            removed["reason"].append("medication-post measure time")
            continue
        
        post_time = post_measure.event_dt.iloc[0] # first post-measure

        if n_meas == "1":
            # first occuring post-measure
            post_times = [post_time]
        elif n_meas == "m":
            # check all post_times until a drug change occurs if multiple measures
            post_times = list(post_measure.event_dt)

        # QC4 & QC5: medication changes & prescription regularity until after post measure
        qc, compl, dose_avg, reason, max_post_time = med_qc(medication_person_df, start_time, drug_start_info["drug"], drug_start_info["dose"], post_times, filter)

        # mitigate the case where missingness is too high because of extended post-time if multiple measures
        if (reason == "prescription missingness too high") & (n_meas == "m"):
            # redo with the first post_time
            post_times = [post_measure.event_dt.iloc[0]]
            qc, compl, dose_avg, reason, max_post_time = med_qc(medication_person_df, start_time, drug_start_info["drug"], drug_start_info["dose"], post_times, filter)

        if qc == False:
            removed["eid"].append(eid)
            removed["reason"].append(reason)
            continue

        post_measures = post_measure[post_measure.event_dt <= max_post_time]
        
        result["eid"].append(eid)
        result["drug_start_date"].append(start_time)
        result["baseline_measure"].append(np.mean(baseline_measure[m])) # average baseline measure
        result["baseline_DT"].append((start_time - baseline_time).days)
        result["post_measure"].append(np.mean(post_measures[m])) # average post measures
        result["post_measure_DT"].append((post_time - start_time).days)
        result["drug_start"].append(drug_start_info["drug"])
        result["dose_start"].append(drug_start_info["dose"])
        result["dose_avg"].append(dose_avg)
        result["completeness"].append(compl)
        result["antilipemic_notstatin_prior"].append(antilipemic_notstatin_prior)

    result_df = pd.DataFrame(result)
    removed_df = pd.DataFrame(removed)

    return result_df, removed_df

#### determine PGx phenotype
result_df, removed_df = PGx_phenotype(medication_df, measure_df, primary_data_df, filter, n_meas, baseline_dt_before, baseline_dt_after, post_dt_min, post_dt_max, prim_care_time)

#### add covariates
cov_df = pd.read_csv(snakemake.input["cov"], sep = '\t')
result_df = result_df.merge(cov_df, on = "eid")
result_df["day"] = 1 # dummy birth day
result_df = result_df.rename(columns = {"birth_year": "year", "birth_month": "month"})
result_df["birth_date"] = pd.to_datetime(result_df[["year", "month", "day"]])
result_df["age_drug"] = ((pd.to_datetime(result_df["drug_start_date"]) - pd.to_datetime(result_df["birth_date"])).dt.days/365.25).astype(int)

result_df.to_csv(snakemake.output[0], sep = '\t', index = False)
removed_df.to_csv(snakemake.output[1], sep = '\t', index = False)