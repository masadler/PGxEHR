library(dplyr)
library(stringr)
library(tidyr)
library(data.table)

m = snakemake@wildcards[["measure"]]

print("Loading clinical data")
clinical_data <- read.table(snakemake@input[["gp_clinical"]], header=T, stringsAsFactors = F, sep="\t", quote="", comment.char = "")
print("clinical data loaded")

code_df <- read.csv(snakemake@input[["codes"]], sep = '\t')
code_df <- code_df[!duplicated(code_df$read_code), c("read_code", "read_term")]

# select measures for read2 data
read2_data <- clinical_data %>% filter(read_2 != "")
read2_data_m <- merge(read2_data, code_df, by.x = "read_2", by.y = "read_code")
read2_data_m <- rename(read2_data_m, read_code = read_2) %>%
                  select(-read_3)

# select measures for ctv3
read3_data <- clinical_data %>% filter(read_3 != "")
read3_data_m <- merge(read3_data, code_df, by.x = "read_3", by.y = "read_code")
read3_data_m <- rename(read3_data_m, read_code = read_3)  %>%
                  select(-read_2)

# extract values
data_m <- rbind(read2_data_m, read3_data_m)
data_m[m] <- as.numeric(data_m$value1)
data_m[(data_m$data_provider == 2) & (data_m$value2 != ""), m] <-  as.numeric(data_m[(data_m$data_provider == 2) & (data_m$value2 != ""), "value2"]) # for data provider 2 -> Scotland

# for SBP and certain codes, take second value (first one being DBP)
if (m == "SBP"){
    codes = c("2464.", "2466.", "2465.", "246d.", "2467.", "246K.")
    data_m[data_m$read_2 %in% codes, m] = as.numeric(data_m[data_m$read_2 %in% codes, "value2"])
    data_m[data_m$read_3 %in% codes, m] = as.numeric(data_m[data_m$read_3 %in% codes, "value2"])
}

data_m <- data_m[!is.na(data_m[, m]), ]

codes <- unique(data_m$read_code)

for (code in codes){
    print(code)
    data_m_code <- data_m[data_m$read_code == code, ]
    print(dim(data_m_code))
    print(summary(as.numeric(data_m_code[, m])))
}

data_m <- data_m[, c("eid", "event_dt", "read_code", "read_term", m)]

if (m == "HbA1c"){
    data_m <- data_m[data_m$read_code %in% c("42W4.", "42W5.", "XaERp", "XaPbt"), ]
    data_m <- data_m[((data_m$read_code %in% c("42W4.", "XaERp")) &
                    (data_m[m] >= 3.9) & (data_m[m] <= 20)) | 
                    ((data_m$read_code %in% c("42W5.", "XaPbt")) &
                    (data_m[m] >= 20) & (data_m[m] <= 195)), ]

    # convert percentage to mmol/mol: (A1c_perc-2.152)/0.09148
    data_m[data_m$read_code %in% c("42W4.", "XaERp"), m] = (data_m[data_m$read_code %in% c("42W4.", "XaERp"), m] - 2.152)/0.09148
}

# transform date
data_m$event_dt = as.Date(data_m$event_dt, '%d/%m/%Y')

# add data from UKBB
header <- as.data.frame(fread(snakemake@input[["pheno"]], header = T, nrow = 0))
pheno_id <- snakemake@params[["m"]]
# get multiple instances
instances <- c("0", "1")
if (m %in% c("SBP", "HR")){
    instances = c(instances, "2", "3")
}

for (inst in instances){
    pheno_df <- as.data.frame(fread(snakemake@input[["pheno"]], header = T, select = c("eid", paste0(pheno_id, "-", inst, ".0")), colClass = "character"))
    names(pheno_df) <- c("eid", m)
    pheno_df <- pheno_df[pheno_df[m] != "", ]

    date_df <- as.data.frame(fread(snakemake@input[["inst"]], header = T, select = c("eid", paste0("53-", inst, ".0")), colClass = "character"))
    names(date_df) <- c("eid", "event_dt")
    date_df$event_dt <- as.Date(date_df$event_dt, '%Y-%m-%d')
    pheno_df <- merge(pheno_df, date_df, by = "eid")

    pheno_df$read_code <- "UKBB"
    pheno_df$read_term <- "UKBB"

    data_m <- rbind(data_m, pheno_df[, c("eid", "event_dt", "read_code", "read_term", m)])
}

write.table(data_m, snakemake@output[[1]], sep = '\t', quote = F, row.names = F)