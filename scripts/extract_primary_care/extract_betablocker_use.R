library(optparse)
library(dplyr)
library(stringr)
library(tidyr)

read_2_lkps <- read.table(snakemake@input[["read_v2"]], stringsAsFactors = F, header=T, sep="\t", quote='', comment.char = '')
names(read_2_lkps) <- c("read_2", "term_description", "status_flag")
read_2_lkps$term_description <- tolower(read_2_lkps$term_description)

print("Loading drug data")
drug_data <- read.table(snakemake@input[["gp_scripts"]], header=T, stringsAsFactors = F, sep="\t", quote="", comment.char = "")
print("Drug data loaded")

drug_data$drug_name <- tolower(enc2utf8(as.character(drug_data$drug_name)))
drug_data$quantity <- tolower(enc2utf8(as.character(drug_data$quantity)))

# eliminate dots in BNF code
drug_data$bnf_code <- gsub("\\.", "", drug_data$bnf_code)

# read drug vocabularies
# data_provider 1 (England Vision) -> read_v2 & dmd
# data_provider 2 (Scotland) -> read_v2 & bnf
# data_provider 3 (England TPP) -> bnf (with dots)
# data_provider 4 (Wales) -> read_v2 (no name & quantity)

#### select betablocker entries ####

print("Extract betablocker entries")

# non-selective: nadolol, pindolol, propranolol, sotalol, timolol
# selective: acebutolol, atenolol, metoprolol, betaxolol, bisoprolol, esmolol, nebivolol
# with alpha: carvedilol, labetolol

# read_2: based on code
bb_read2_df = read.table(snakemake@input[["read2"]], header = T)
#bb_read2_df = read.table("data/betablocker_classification_read2.tsv", header = T)
bb_read2_df$category = "betablocker"

drug_data_read2 <- drug_data %>% 
                  filter(read_2 != "") %>%
                  mutate(code = str_sub(read_2, 1, 3))
drug_data_read2 <- merge(drug_data_read2, bb_read2_df)

# If drug name for read_2 is missing -> replace with description from read_2_lkp
read_code_data <- drug_data_read2 %>% 
  select(read_2) %>% 
  unique %>% 
  arrange %>% 
  separate(read_2, sep="\\.", into=c('read_2', 'trailing')) %>% 
  mutate(read_2 = paste(read_2, '.', sep='')) %>%
  select(-trailing) %>%
  unique

read_code_data <- left_join(read_code_data, read_2_lkps) %>% 
  filter(!is.na(term_description)) %>% 
  select(-status_flag)

# fix the read code in the drug data
drug_data_read2 <- drug_data_read2 %>%
  separate(read_2, sep="\\.", into=c('read_2', 'trailing')) %>% 
  mutate(read_2 = paste(read_2, '.', sep='')) %>%
  select(-trailing)

drug_data_read2 <- left_join(drug_data_read2, read_code_data)

drug_data_read2 <- drug_data_read2 %>% 
  mutate(drug_name = case_when(
    !is.na(term_description) ~ term_description,
    is.na(term_description) ~ drug_name))

# Extract beta-blockers from BNF
# pre-filter on Beta-Adrenoceptor Blocking Drugs
drug_data_bnf <- drug_data %>% filter(str_starts(bnf_code, "02040"))
drug_data_bnf$category <- "betablocker"
drug_data_bnf$drug <- as.character(NA)

# do a string search per beta blocker drug # 02.04.0
bb_bnf_df = read.csv(snakemake@input[["bnf"]], sep = '\t')
#bb_bnf_df = read.csv("data/betablocker_bnf_keyword_search.tsv", sep = '\t')
bb_bnf_df$keyword = tolower(bb_bnf_df$keyword)
drugs = as.vector(unique(bb_bnf_df$drug))

for (d in drugs){
  code = "02040"
  kws = paste(as.vector(bb_bnf_df[bb_bnf_df$drug == d, "keyword"]),collapse = '|')

  drug_data_bnf <- drug_data_bnf %>% 
                      mutate(drug = case_when(
                        !is.na(drug) ~ drug,
                        (is.na(drug) & str_detect(drug_name, kws) & str_starts(bnf_code, code)) ~ d
                        )) 
}

# concatenate data and select columns of interest
cols = c("eid", "issue_date", "drug_name", "quantity", "drug", "category")
drug_data <- rbind(drug_data_read2[, cols], drug_data_bnf[, cols])

# mark beta-blocker combination therapies 
combis = c("clopamide", "viskaldix", "hydchloroth", "sotazide", "tolerzide",
                "gppe", "bendroflumeth", "moducren", "prestim", "monozide", "tenben",
                "aspirin", "secadrex", "kalten", "spironol", "inderetic", "inderex", 
                "spiroprop", "co-betaloc", "lopresoretic", "corgaretic")
drug_data$combi <- as.character(NA)

for (i in 1:length(combis)){
  c = combis[i]
  
  drug_data = drug_data %>% 
                    mutate(combi = case_when(
                      ((str_detect(drug_name, c)) & (category == "betablocker")) ~ c,
                       !is.na(combi) ~ combi))
  i = i+1
}

#### manual corrections

## combinations with beta-blocker atenolol
drug_combis = c("beta-adalat", "co-tenidone", "kalten", "tenif", "tenoret", "atenolol+nifedipine", "atenolol+bendrofluazide",
                "atenolol+co-amilozide")

for (drug_combi in drug_combis){
  drug_data = drug_data %>% 
                mutate(drug = ifelse(str_detect(drug_name, drug_combi), "atenolol", as.character(drug)),
                        category = ifelse(str_detect(drug_name, drug_combi), "betablocker", as.character(category)),
                        combi = ifelse(str_detect(drug_name, drug_combi), drug_combi, as.character(combi)))

}

## combinations with beta-blocker metoprolol
drug_combis = c("metoprolol tart+hydrochlorothiazide", "propranolol hcl+bendrofluazide", "sotalol hcl+hydrochlorothiazide")
drugs = c("metoprolol", "propranolol", "sotalol")

for (i in 1:length(drug_combis)){
  drug_combi = drug_combis[i]
  drug_ = drugs[i]
  drug_data = drug_data %>% 
                mutate(drug = ifelse(str_detect(drug_name, drug_combi), drug_, as.character(drug)),
                        category = ifelse(str_detect(drug_name, drug_combi), "betablocker", as.character(category)),
                        combi = ifelse(str_detect(drug_name, drug_combi), drug_combi, as.character(combi)))

}

####  extract dose and quantity ####

# Function to determine the dosage for each row in the script data
determine_dose <- function(drug_data) {
  
  # calculate dose by calculating total mg in prescription
  print("Calculating dose")
  drug_data <- drug_data %>%
    mutate(drug_name = tolower(drug_name)) %>%
    mutate(dose = case_when(
      str_detect(drug_name, "\\d+.?\\d+\\s?mg") ~ str_extract(drug_name, "\\d+.?\\d+\\s?mg"),
      str_detect(drug_name, "\\d+\\s?mg") ~ str_extract(drug_name, "\\d+\\s?mg"),
      str_detect(drug_name, "\\d+.?\\d+\\s?microgram") ~ str_extract(drug_name, "\\d+.?\\d+\\s?microgram"),
      str_detect(drug_name, "\\d+\\s?microgram") ~ str_extract(drug_name, "\\d+\\s?microgram"),
      str_detect(drug_name, "\\d+.?\\d+\\s?mcg") ~ str_extract(drug_name, "\\d+.?\\d+\\s?mcg"),
      str_detect(drug_name, "\\d+\\s?mcg") ~ str_extract(drug_name, "\\d+\\s?mcg"),
      str_detect(drug_name, "\\d+.?\\d+\\s?nanogram") ~ str_extract(drug_name, "\\d+.?\\d+\\s?nanogram"),
      str_detect(drug_name, "\\d+\\s?nanogram") ~ str_extract(drug_name, "\\d+\\s?nanogram"),
      str_detect(drug_name, "\\d+.?\\d+\\s?milligram") ~ str_extract(drug_name, "\\d+.?\\d+\\s?milligram"),
      str_detect(drug_name, "\\d+\\s?milligram") ~ str_extract(drug_name, "\\d+\\s?milligram"),
      str_detect(drug_name, "\\d+.?\\d+\\s?g") ~ str_extract(drug_name, "\\d+.?\\d+\\s?g"),
      str_detect(drug_name, "\\d+\\s?g") ~ str_extract(drug_name, "\\d+\\s?g")
    )) %>% 
    mutate(dose = case_when(
      str_detect(dose, "microgram") ~ (as.numeric(str_extract(drug_name, "\\d+")) / 1000),
      str_detect(dose, "mcg") ~ (as.numeric(str_extract(drug_name, "\\d+")) / 1000),
      str_detect(dose, "nanogram") ~ (as.numeric(str_extract(drug_name, "\\d+")) / 1000000),
      str_detect(dose, "mg") ~ (as.numeric(str_extract(drug_name, "[[:digit:]]+\\.*[[:digit:]]*"))), 
      str_detect(dose, "g") ~ (as.numeric(str_extract(drug_name, "\\d+")) * 1000), # metformin has grams
    )) 
  
 # Add quantity
  print("Calculating quantity")
  drug_data <- drug_data %>%
    mutate(packs = tolower(quantity)) %>%
    mutate(packs = case_when(str_detect(quantity, 'pack') ~ str_extract(quantity, "\\d+\\s?pack"))) %>%
    mutate(packs = str_replace(packs, 'pack', "")) %>%
    mutate(quantity = str_replace(quantity, '\\(', "")) %>%
    mutate(quantity = str_replace(quantity, '\\)', "")) %>%
    mutate(quantity = str_replace(quantity, '\\[', "")) %>%
    mutate(quantity = str_replace(quantity, '\\]', "")) %>%
    mutate(quantity = str_replace(quantity, 'x', "")) %>%
    mutate(amount = case_when(
      str_detect(quantity, "cap") ~ str_extract(quantity, "\\d+\\s?-?\\s?cap"),
      str_detect(quantity, "boot") ~ str_extract(quantity, "\\d+\\s?-?\\s?boot"),
      str_detect(quantity, "cp") ~ str_extract(quantity, "\\d+\\s?-?\\s?cp"),
      str_detect(quantity, "tab") ~ str_extract(quantity, "\\d+\\s?-?\\s?tab"),
      str_detect(quantity, "x tablets") ~ str_extract(quantity, "\\d+\\s?x"),
      str_detect(quantity, "pack of") ~ str_extract(quantity, "pack of \\d+"),
      str_detect(quantity, "unit") ~ str_extract(quantity, "\\d+\\s?unit"),
      str_detect(quantity, "day") ~ str_extract(quantity, "\\d+\\s?day"),
      str_detect(quantity, "sachet") ~ str_extract(quantity, "\\d+\\s?sach"),
      str_detect(quantity, "millilitres") ~ str_extract(quantity, "\\d+\\s?millilitres"),
      str_detect(quantity, "ml") ~ str_extract(quantity, "\\d+\\s?ml"),
      T ~ quantity
    )) %>%
    mutate(amount = str_replace_all(amount, fixed(" "), "")) %>%
    mutate(amount = str_replace(amount, '\\(', "")) %>%
    mutate(amount = str_replace(amount, '\\)', "")) %>%
    mutate(amount = str_replace(amount, 'cap', "")) %>%
    mutate(amount = str_replace(amount, 'boot', "")) %>%
    mutate(amount = str_replace(amount, '-', "")) %>%
    mutate(amount = str_replace(amount, 'cp', "")) %>%
    mutate(amount = str_replace(amount, 'tab', "")) %>%
    mutate(amount = str_replace(amount, 'x', "")) %>%
    mutate(amount = str_replace(amount, 'ml', "")) %>%
    mutate(amount = str_replace(amount, 'millilitres', "")) %>%
    mutate(amount = str_replace(amount, 'day', "")) %>%
    mutate(amount = str_replace(amount, 'unit', "")) %>%
    mutate(amount = str_replace(amount, 'packof', "")) %>%
    mutate(amount = str_replace(amount, 'sach', "")) %>%
    mutate(amount = as.numeric(amount), packs=as.numeric(packs), dose=as.numeric(dose)) %>%
    # misc extra cases i've seen
    mutate(amount=case_when(
      quantity == '1 150' ~ 150,
      quantity == '60 - 100 mg' ~ 60,
      quantity == '1 x 100' ~ 1,
      quantity == '112 - 100 mg' ~ 112,
      amount < 0 ~ 0,
      str_detect(quantity, "[*]28 tab") ~ as.numeric(str_extract(quantity, "\\d+"))*as.numeric(str_extract(str_extract(quantity, "[*]\\d+"), "\\d+")),
      quantity == '4  28 tablet' ~ 112,
      T~amount
    )) %>%
    mutate(packs = case_when(
      is.na(packs)~1,
      T~packs
    ))

  # detect volume measurements
  drug_data <- drug_data %>%
    mutate(volume = str_detect(drug_name, "\\d+.?\\d+\\s?mg\\s?/\\s?\\d+\\s?ml")) %>%
    mutate(fraction = str_detect(drug_name, "\\d+.?\\d+\\s?%"))
  #  filter(fraction == F, volume == F)
    
  file=snakemake@output[["no_amount"]]
  no_amount <- drug_data[is.na(drug_data$amount),] %>% group_by(quantity) %>% tally %>% arrange(-n)
  write.table(no_amount, file, row.names = F, sep="\t", quote=F)      
   
  file=snakemake@output[["no_dose"]]
  no_dose <- drug_data[is.na(drug_data$dose),] %>% group_by(drug_name) %>% tally %>% arrange(-n)
  write.table(no_dose, file, row.names = F, sep="\t", quote=F)    

  drug_data <- drug_data %>%
    #filter(!is.na(packs), !is.na(amount), !is.na(dose)) %>%
    #filter(!is.na(dose)) %>% # allow to have NAs on amount & packs
    mutate(prescription_mg = packs*amount*dose)
  
  return(drug_data)
}

#################################################

# format the date, sort, remove any without a date, remove duplicates of data & drug
# removing duplicates will remove some real multi-prescriptions, but that would not be
# relevant here as it would be the exact same drug, dose, etc.
drug_data <- drug_data %>% 
  mutate(issue_date_new = as.Date(issue_date, '%d/%m/%Y')) %>%
  arrange(issue_date_new) %>%
  filter(!is.na(issue_date_new)) %>%
  mutate(issue_date=issue_date_new) %>%
  select(-issue_date_new) %>%
  arrange(eid, issue_date) %>%
  unique

drug_data <- determine_dose(drug_data)

file=snakemake@output[["dose"]]
write.table(drug_data, file, row.names = F, sep="\t", quote=F)