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

#### select metformin entries ####

print("Extract metformin entries")

# read_2: based on code

# metformin
drug_data_read2 <- filter(drug_data, grepl("^f4|^f3|^f2|^f1|^ft", read_2))

drug_data_read2$category <- as.character(NA)

cats <- c("biguanides", "sulfonylureas", "insulin_short", "insulin", "other")
codes <- c("f4", "f3", "f1", "f2", "ft") # all f4 are metformin (f41)

i = 1
for (cat in cats){

  drug_data_read2 = drug_data_read2 %>% 
                    mutate(category = case_when(
                      (str_starts(read_2, codes[i])) ~ cat,
                       !is.na(category) ~ category)) 
  i = i+1
}

drug_data_read2$drug <- as.character(NA)
drug_data_read2[(drug_data_read2$cat == "biguanides") & (!is.na(drug_data_read2$cat)), "drug"] = "metformin"

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

# Extract metformin from BNF

# pre-filter for antidiabetic
drug_data_bnf <- filter(drug_data, grepl("^060102|^060101", bnf_code))
drug_data_bnf$drug <- as.character(NA)
drug_data_bnf$category <- as.character(NA)

# annotate with drug and category

# "06.01.02.02.00" -> metformin & others
# "06.01.02.01.00" -> sulfonylureas
# "06.01.02.03.00" -> others
# "06.01.01" -> insulin

cats <- c("biguanides", "sulfonylureas", "other", "insulin_pen", "insulin_short", "insulin", "insulin_pen")
codes <- c("0601020200", "0601020100", "0601020300", "0601010300", "0601010100",
           "0601010200", "0601010400")

combi = paste(c("sitagliptin", "janumet", "linagliptin", "jentadueto", "saxagliptin",
                "komboglyze", "alogliptin", "vipdomet", "dapagliflozin", "xigduo", "canagliflozin",
                "vokanamet", "empagliflozin", "synjardy", "rosiglitazone",
                "avandamet", "pioglitazone", "competact", "vildagliptin",
                "eucreas"), collapse = '|')

i = 1
for (cat in cats){
  drug_data_bnf[drug_data_bnf$bnf_code == codes[i], "category"] = cat

  if (cat == "biguanides"){
    drug_data_bnf = drug_data_bnf %>% 
                    mutate(category = case_when(
                      ((!is.na(category) & str_detect(drug_name, combi) & (category == cat)) ~ "other"),
                      (!is.na(category) & !str_detect(drug_name, combi) & (category == cat)) ~ cat))
  }
  i = i+1
}

drug_data_bnf[(drug_data_bnf$category == "biguanides") & !is.na(drug_data_bnf$category), "drug"] = "metformin"

# concatenate data and select columns of interest
cols = c("eid", "issue_date", "drug_name", "quantity", "drug", "category")
drug_data <- rbind(drug_data_read2[, cols], drug_data_bnf[, cols])

# annotate combination therapies with metformin
met_combis = paste0(c("metformin", "avandamet", "competact", "eucreas", "janumet", "jentadueto", "komboglyze",
               "vipdomet", "xigduo", "vokanamet", "synjardy"), collapse = '|')

drug_data$combi <- "no"
drug_data[!is.na(drug_data$category) & (drug_data$category == "other") & (str_detect(drug_data$drug_name, met_combis)), "combi"] <- "yes"

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

