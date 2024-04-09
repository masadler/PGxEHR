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

#### select statin entries ####

print("Extract statin entries")

# read_2: based on code

# {"atorvastatin": [1545958], "fluvastatin": [1549686], 
# "lovastatin": [1592085], "pitavastatin": [40165636], 
# "pravastatin": [1551860], "rosuvastatin": [1510813],
# "simvastatin": [1539403]}

# Notes: no lovastatin, pitavastatin in read2; cerivastatin in lkp, but removed from market

antilipemic_read2 <- c("bx") # define a broad antilipemic medication category
drug_data_read2 <- drug_data %>% 
                  filter(read_2 != "") %>%
                  mutate(code2 = str_sub(read_2, 1, 2))
drug_data_read2 <- drug_data_read2[drug_data_read2$code2 %in% antilipemic_read2, ]
drug_data_read2$antilipemic <- "yes"
drug_data_read2$drug <- as.character(NA)

statins <- c("atorvastatin", "fluvastatin", "pravastatin", "rosuvastatin", "simvastatin", "cerivastatin")
codes <- c("bxi", "bxg", "bxe", "bxk", "bxd", "bxj")

i = 1
for (statin in statins){

  drug_data_read2 = drug_data_read2 %>% 
                    mutate(drug = case_when(
                      (str_starts(read_2, codes[i])) ~ statin,
                       !is.na(drug) ~ drug))
  i = i+1
}

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

# Extract statins from BNF

# get all antilipemic
drug_data_bnf <- drug_data %>% filter(str_starts(bnf_code, "02120"))
drug_data_bnf$antilipemic <- "yes"
drug_data_bnf$drug <- as.character(NA)

# do a string search per statin type
statin_bnf_df = read.table(snakemake@input[["bnf"]], header = T)
statin_bnf_df$Keywords = tolower(statin_bnf_df$Keywords)
statins = as.character(unique(statin_bnf_df$Type))

for (i in 1:length(statins)){
  statin = statins[i]
  kws = paste(as.vector(statin_bnf_df[statin_bnf_df$Type == statin, "Keywords"]),collapse = '|')

  drug_data_bnf = drug_data_bnf %>% 
                    mutate(drug = case_when(
                      ((str_detect(drug_name, kws)) & (bnf_code != "")) ~ statin,
                       !is.na(drug) ~ drug))
  i = i+1
}

# concatenate data and select columns of interest
cols = c("eid", "issue_date", "drug_name", "quantity", "drug", "antilipemic")
drug_data <- rbind(drug_data_read2[, cols], drug_data_bnf[, cols])

# mark combination therapies ("colib": fenofibrate/simvastatin , "inegy": "Simvastatin & ezetimibe")
#combis = paste(c("ezetimibe", "inegy", "fenofibrate", "cholib"), collapse = '|')
drug_data$combi <- as.character(NA)
combis = c("inegy", "cholib") # combis that contain statins

for (i in 1:length(combis)){
  c = combis[i]
  
  drug_data = drug_data %>% 
                    mutate(combi = case_when(
                      ((str_detect(drug_name, c))) ~ c,
                       !is.na(combi) ~ combi))
  i = i+1
}

combis = c("ezetimibe", "fenofibrate") # exist with statins in combination

for (i in 1:length(combis)){
  c = combis[i]
  
  drug_data = drug_data %>% 
                    mutate(combi = case_when(
                      ((str_detect(drug_name, c) & str_detect(drug_name, "statin"))) ~ c,
                       !is.na(combi) ~ combi))
  i = i+1
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


# this is not used in the end
determine_statin_strength <- function(statin, dose){
  # assign statin strength: https://compendiumapp.com/user_uploads/000/001/898_BGqLwz_Capture.PNG
  # issue: cerivastatin is missing
    if (is.na(statin)){
      return (c(NA, NA, 0))
    } else if (statin == "fluvastatin"){
        if (is.na(dose)){
            return (c(18, 34, 0))
        } else {
            if (dose <= 20){
                return (c(18, 20, 1))
            } else if ((dose > 20) & (dose <= 40)){
                return (c(20, 24, 1))
            } else if (dose > 40){
                return (c(24, 34, 1))
              }
          }
     } else if (statin == "pravastatin"){
        if (is.na(dose)){
            return (c(18, 36, 0))
        } else {
            if (dose <= 10){
                return (c(18, 20, 1))
            } else if ((dose > 10) & (dose <= 20)){
                return (c(20, 30, 1))
            } else if ((dose > 20) & (dose <= 40)){
               return (c(30, 32, 1))
            } else if (dose > 40){
                return (c(32, 36, 1))
            }
           }
     } else if (statin == "lovastatin"){
        if (is.na(dose)){
            return (c(18, 38, 0))
        } else {
            if (dose <= 10){
                return (c(18, 20, 1))
            } else if ((dose > 10) & (dose <= 20)){
                return (c(20, 26, 1))
            } else if ((dose > 20) & (dose <= 40)){
                return (c(26, 30, 1))
            } else if ((dose > 40) & (dose <= 80)){
                return (c(30, 38, 1))
            }
         }
     } else if (statin == "simvastatin"){
        if (is.na(dose)){
            return (c(18, 40, 0))
        } else {
            if (dose <= 5){
                return (c(18, 26, 1))
            } else if ((dose > 5) & (dose <= 10)){
                return (c(26, 30, 1))
            } else if ((dose > 10) & (dose <= 20)){
                return (c(30, 34, 1))
            } else if (dose > 20){
                return (c(34, 40, 1))
            }
        }
     } else if (statin == "pitavastatin"){
        if (is.na(dose)){
            return (c(18, 42, 0))
        } else {
            if (dose <= 1){
                return (c(18, 30, 1))
            } else if ((dose > 1) & (dose <= 2)){
                return (c(30, 32, 1))
            } else if (dose > 2) {
                return (c(32, 42, 1))
            }
          }
     } else if (statin == "atorvastatin"){
        if (is.na(dose)){
            return (c(18, 60, 0))
        } else {
            if (dose <= 10){
                return (c(18, 38, 1))
            } else if ((dose > 10) & (dose <= 20)){
                return (c(38, 42, 1))
            } else if ((dose > 20) & (dose <= 40)){
                return (c(42, 50, 1))
            } else if ((dose > 40)){
                return (c(50, 60, 1))
            }
        }
     } else if (statin == "rosuvastatin"){
        if (is.na(dose)){
            return (c(18, 63, 0))
        } else {
            if (dose <= 5){
                return (c(18, 44, 1))
            } else if ((dose > 5) & (dose <= 10)){
                return (c(44, 52, 1))
            } else if ((dose > 10) & (dose <= 20)){
                return (c(46, 52, 1))
            } else if (dose > 20){
                return (c(52, 63, 1))
            }
        }
     } else {
        return (c(NA, NA, 0))
     }
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

#### determine statin strength - MIN & MAX #### 

results <- apply(drug_data, 1, function(x) determine_statin_strength(x['drug'], x['dose']))

drug_data$strength_min <- results[1,]
drug_data$strength_max <- results[2,]
drug_data$dose_available <- results[3,]

file=snakemake@output[["dose"]]
write.table(drug_data, file, row.names = F, sep="\t", quote=F)   



