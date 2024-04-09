library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

plot_cohort <- function(cohort, drug, biomarker, ylab, colors = c("#cfdde6", "#3778a3")){

    df = read.csv(paste0("output/PGx_phenotype/drug_", cohort, "_l_m_PGx_phenotype.tsv"), sep = "\t", header = T)

    if (cohort == "statin_LDL"){
        unit = "mmol/L"
        l = 2
    } else if (cohort == "statin_TC"){
        unit = "mmol/L"
        l = 5
    } else if (cohort == "statin_HDL"){
        unit = "mmol/L"
        l = 0.2
        } else if (cohort == "metformin_HbA1c"){
        unit = "mmol/mol"
        l = 42
    } else if (cohort %in% c("antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
        "antihpt_betablocker_SBP")){
        unit = "mmHg"
        l = 120
    } else if (cohort == "betablocker_HR"){
        unit = "beats/min"
        l = 60
    }

    df = df[df$baseline_measure > l,]

    med_count = df %>% count(drug_start)
    drugs = as.vector(med_count[med_count$n >= 20, "drug_start"])
    df = df[df$drug_start %in% drugs, ]

    df$cohort = drug
    df = df[, c("cohort", "baseline_measure", "post_measure")]
    names(df) = c("cohort", "measure1", "measure2")
    long_df = df %>% pivot_longer(!cohort, names_to = "measure", values_to = "level")

    control_df = read.csv(paste0("output/PGx_phenotype/control_", cohort, "_l_1_PGx_phenotype.tsv"), sep = "\t", header = T)
    control_df$cohort = "control"
    control_df = control_df[, c("cohort", "measure1", "measure2")]
    long_df_control = control_df %>% pivot_longer(!cohort, names_to = "measure", values_to = "level")
    
    df = rbind(long_df, long_df_control)

    df$biomarker = biomarker

    df$biomarker = factor(df$biomarker)
    df$cohort = factor(df$cohort, levels = c(drug, "control"))
    df$measure = factor(df$measure, levels = c("measure1", "measure2"))

    p1 = ggplot(df, aes(x = cohort, y = level)) +
    geom_boxplot(aes(fill = measure)) + 
    facet_wrap(~biomarker) +
    labs(x = "", y = ylab, fill = "") +
    scale_fill_manual(values = colors) +
    theme(legend.position="none")

    return(p1)
}

# single drug-biomarker combination
cohorts = c("statin_LDL", "statin_TC", "statin_HDL", "metformin_HbA1c", "betablocker_HR")
drugs = c("statin", "statin", "statin", "metformin", "beta blocker")
biomarkers = c("LDL-C", "TC", "HDL-C", "HbA1c", "Heart rate")
units = c("mmol/L", "mmol/L", "mmol/L", "mmol/mol", "beats/min")
colors = list(c("#045291", "#72b3e8"), c("#045291", "#72b3e8"), c("#045291", "#72b3e8"),
              c("#993808", "#e38f66"), c("#7d0575", "#e66cdd"))

for (i in 1:length(cohorts)){
    cohort = cohorts[i]
    drug = drugs[i]
    biomarker = biomarkers[i]
    unit = units[i]
    color = colors[[i]]

    pdf(paste0("output/Figures/cohort_", cohort, "_and_control_before_after_boxplot.pdf"), width = 2.5, height = 2.08)
    print(plot_cohort(cohort, drug, biomarker, unit, colors = color))
    dev.off()

}

# all antihypertensives
plot_cohort_SBP <- function(colors = c("#038216", "#89f098")){

    control_df = read.csv(paste0("output/PGx_phenotype/control_antihpt_all_SBP_l_1_PGx_phenotype.tsv"), sep = "\t", header = T)
    control_df$cohort = "control"
    control_df = control_df[, c("cohort", "measure1", "measure2")]
    long_df = control_df %>% pivot_longer(!cohort, names_to = "measure", values_to = "level")

    cohorts = c("antihpt_all_SBP", "antihpt_ACEi_SBP", "antihpt_CCB_SBP", "antihpt_thiazide_diuretics_SBP", 
                "antihpt_betablocker_SBP")
    cohort_names = c("All (first-line)", "ACEi", "CCB", "diuretics", "beta blocker")

    l = 120

    for (i in 1:length(cohorts)){
        cohort = cohorts[i]
        df = read.csv(paste0("output/PGx_phenotype/drug_", cohort, "_l_m_PGx_phenotype.tsv"), sep = "\t", header = T)
        df$cohort = cohort_names[i]
        df = df[df$baseline_measure > l,]

        med_count = df %>% count(drug_start)
        drugs = as.vector(med_count[med_count$n >= 20, "drug_start"])
        df = df[df$drug_start %in% drugs, ]

        df = df[, c("cohort", "baseline_measure", "post_measure")]
        names(df) = c("cohort", "measure1", "measure2")
        long_df_drug = df %>% pivot_longer(!cohort, names_to = "measure", values_to = "level")

        long_df = rbind(long_df, long_df_drug)
    }

    df = long_df
    
    df$biomarker = "Systolic blood pressure"

    df$biomarker = factor(df$biomarker)
    df$cohort = factor(df$cohort, levels = c("All", "ACEi", "CCB", "diuretics", "beta blocker", "control"))
    df$measure = factor(df$measure, levels = c("measure1", "measure2"))

    p1 = ggplot(df, aes(x = cohort, y = level)) +
    geom_boxplot(aes(fill = measure)) + 
    facet_wrap(~biomarker) +
    labs(x = "", y = "mmHg", fill = "") +
    scale_fill_manual(values = colors) +
    theme(legend.position="none")

    return(p1)
}

pdf(paste0("output/Figures/cohort_SBP_and_control_before_after_boxplot.pdf"), width = 7.5, height = 2.08)
print(plot_cohort_SBP())
dev.off()