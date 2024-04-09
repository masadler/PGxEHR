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
# data_provider 3 (England TPP) -> bnf 
# data_provider 4 (Wales) -> read_v2 (no name & quantity)

#### select anti-hypertensive entries ####

print("Extract anti-hypertensive entries")

# read_2: based on code
bb_read2_df = read.csv(snakemake@input[["read2_antih"]], sep = '\t')

# add beta-blocker information
bb_read2_df_2 = read.table(snakemake@input[["read2_bb"]], header = T)
# change type to "betablocker"
bb_read2_df_2$type = "betablocker"

bb_read2_df = rbind(bb_read2_df, bb_read2_df_2)

# other drugs with anti-hypertensive properties
# b2: thiazide diuretics; b3: loop diuretics; b4: potassium sparing diuretics
# b5: potassium sparing compound diuretics; b6: osmotic diuretics; b7: mercurial diuretics
# b9: diuretics + potassium supplement, bA: calcium channel blocker + angiotensin-converting enzyme inhibitor
# bd: beta-adrenoceptor blockers; be: vasodilator antihypertensives; bf: central antihypertensives
# bg: adrenergic neurone blockers; bh: alpha-adrenoceptor blockers; bi: angiotensin-converting enzyme inhibitors
# bj: ganglion blocking drugs; bk: other antihypertensives

antihypertensives_read2 <- c("b2", "b3", "b4", "b5", "b6", "b7", "b9", "bA", "bd", "be", "bf", "bg", "bh", "bi", "bj", "bk")

drug_data_read2 <- drug_data %>% 
                  filter(read_2 != "") %>%
                  mutate(code2 = str_sub(read_2, 1, 2))
drug_data_read2 <- drug_data_read2[drug_data_read2$code2 %in% antihypertensives_read2, ]
drug_data_read2$category <- "antihypertensive"

drug_data_read2 <- drug_data_read2 %>% 
                  filter(read_2 != "") %>%
                  mutate(code = str_sub(read_2, 1, 3))
drug_data_read2 <- merge(drug_data_read2, bb_read2_df, all.x = T)

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

# Extract anti-hypertensives from BNF
# get codes
bnf_lkps_df = read.csv(snakemake@input[["bnf_all"]], sep = '\t')
bnf_lkps_df$BNF_Presentation_Code = as.character(bnf_lkps_df$BNF_Presentation_Code)
bnf_lkps_df = filter(bnf_lkps_df, grepl("^020201|^020202|^020203|^020204|^020205|^020208|^02040|^02050|^020602", BNF_Presentation_Code))
bnf_lkps_df$BNF_Chemical_Substance = tolower(bnf_lkps_df$BNF_Chemical_Substance)
bnf_lkps_df$drug = bnf_lkps_df$BNF_Chemical_Substance

# classes
classes <- c("Vasodilator Antihypertensive Drugs", "Centrally-Acting Antihypertensive Drugs",
             "Adrenergic Neurone Blocking Drugs", "Alpha-Adrenoceptor Blocking Drugs", "Angiotensin-Converting Enzyme Inhibitors",
             "Angiotensin-II Receptor Antagonists", "Renin Inhibitors", "Calcium-Channel Blockers",
             "Loop Diuretics", "Pot-Sparing Diuretics&Aldosterone Antag")

antihyper_bnf_df = bnf_lkps_df[bnf_lkps_df$BNF_Subparagraph == "Thiazides And Related Diuretics", ]
antihyper_bnf_df$type = "thiazide_diuretics"

for (class in classes){
  print(paste0("\n", class))
  bnf_lkps_df_class = bnf_lkps_df[bnf_lkps_df$BNF_Subparagraph == class, ]

  if (class == "Vasodilator Antihypertensive Drugs"){
    drugs = c("riociguat", "macitentan", "diazoxide", "hydralazine hydrochloride", "minoxidil", "sodium nitroprusside", "bosentan",
    "iloprost", "sitaxentan sodium", "ambrisentan", "sildenafil(vasodilator antihypertensive)",
    "tadalafil (vasodilator antihypertensive)")
    bnf_lkps_df_class = bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance %in% drugs, ]
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "sildenafil(vasodilator antihypertensive)", "drug"] = "sildenafil"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "tadalafil (vasodilator antihypertensive)", "drug"] = "tadalafil"
    bnf_lkps_df_class$type = "vasodilator"
  } else if (class == "Centrally-Acting Antihypertensive Drugs") {
    drugs = c("clonidine hydrochloride", "guanfacine hydrochloride", "methyldopa", "methyldopate hydrochloride", "moxonidine")
    bnf_lkps_df_class = bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance %in% drugs, ]
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "methyldopate", "drug"] = "methyldopa"
    bnf_lkps_df_class$type = "alpha_agonist"
  } else if (class == "Adrenergic Neurone Blocking Drugs") {
    bnf_lkps_df_class$type = "adrenergic_neurone_blocker"
  } else if (class == "Alpha-Adrenoceptor Blocking Drugs") {
    drugs = c("doxazosin mesilate", "indoramin", "phenoxybenzamine hydrochloride", "phentolamine mesilate", 
              "prazosin hydrochloride", "terazosin hydrochloride")
    bnf_lkps_df_class = bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance %in% drugs, ]
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "doxazosin mesilate", "drug"] = "doxazosin"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "phenoxybenzamine hydrochloride", "drug"] = "phenoxybenzamine"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "phentolamine mesilate", "drug"] = "phentolamine"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "prazosin hydrochloride", "drug"] = "prazosin"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "terazosin hydrochloride", "drug"] = "terazosin"
    bnf_lkps_df_class$type = "alpha_blocker"
  } else if (class == "Angiotensin-Converting Enzyme Inhibitors") {
    drugs = c("perindopril tosilate", "moexipril hydrochloride", "cilazapril", "captopril", 
              "enalapril maleate", "fosinopril sodium", "lisinopril", "perindopril erbumine",
              "quinapril hydrochloride", "ramipril", "trandolapril", "imidapril hydrochloride",
              "benazepril hydrochloride", "perindopril arginine")
    bnf_lkps_df_class = bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance %in% drugs, ]
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "perindopril tosilate", "drug"] = "perindopril"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "moexipril hydrochloride", "drug"] = "moexipril"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "enalapril maleate", "drug"] = "enalapril"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "fosinopril sodium", "drug"] = "fosinopril"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "quinapril hydrochloride", "drug"] = "quinapril"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "imidapril hydrochloride", "drug"] = "imidapril"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "benazepril hydrochloride", "drug"] = "benazepril"
    bnf_lkps_df_class$type = "ACEi"
  } else if (class == "Angiotensin-II Receptor Antagonists") {
    drugs = c("azilsartan medoxomil", "olmesartan medoxomil", "candesartan cilexetil", "irbesartan", 
              "losartan potassium", "telmisartan", "valsartan", "eprosartan")
    bnf_lkps_df_class = bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance %in% drugs, ]
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "azilsartan medoxomil", "drug"] = "azilsartan"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "olmesartan medoxomil", "drug"] = "olmesartan"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "candesartan cilexetil", "drug"] = "candesartan"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "losartan potassium", "drug"] = "losartan"
    bnf_lkps_df_class$type = "ARB"
  } else if (class == "Renin Inhibitors") {
    bnf_lkps_df_class$type = "renin-i"
  } else if (class == "Calcium-Channel Blockers") {
    notdrugs = c("diltiazem hydrochloride with thiazides", "valsartan/amlodipine", "buflomedil hydrochloride", 
                  "perhexiline maleate", "trimetazidine hydrochloride", "prenylaminelactate")
    bnf_lkps_df_class = bnf_lkps_df_class[!(bnf_lkps_df_class$BNF_Chemical_Substance %in% notdrugs), ]
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "diltiazem hydrochloride", "drug"] = "diltiazem"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "lercanidipine hydrochloride", "drug"] = "lercanidipine"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "nicardipine hydrochloride", "drug"] = "nicardipine"
    bnf_lkps_df_class[bnf_lkps_df_class$BNF_Chemical_Substance == "verapamil hydrochloride", "drug"] = "verapamil"
    bnf_lkps_df_class$type = "CCB"
  } else if (class == "Loop Diuretics"){
    notdrugs = c("other loop diuretic preps")
    bnf_lkps_df_class = bnf_lkps_df_class[!(bnf_lkps_df_class$BNF_Chemical_Substance %in% notdrugs), ]
    bnf_lkps_df_class$type = "loop_diuretics"
  } else if (class == "Pot-Sparing Diuretics&Aldosterone Antag"){
    bnf_lkps_df_class$type = "pot_sparing_diuretics"
  }
  antihyper_bnf_df = rbind(antihyper_bnf_df, bnf_lkps_df_class)
}

antihyper_bnf_df <- antihyper_bnf_df %>%
                    mutate(root = tolower(str_split(BNF_Presentation, "_| ", simplify = TRUE)[, 1]))

# add code
antihyper_bnf_code_df <- data.frame(code = c("020201", "020501", "020502", 
                                             "020503", "020504", "02050501",
                                            "02050502", "02050503", "020602", 
                                            "020202", "020203"), 
                                   type = c("thiazide_diuretics", "vasodilator", "alpha_agonist", 
                                            "adrenergic_neurone_blocker", "alpha_blocker", "ACEi", 
                                            "ARB", "renin-i", "CCB", 
                                            "loop_diuretics", "pot_sparing_diuretics"))
antihyper_bnf_df <- merge(antihyper_bnf_df, antihyper_bnf_code_df, by = "type")

# pre-filter for anti-hypertensives
drug_data_bnf <- filter(drug_data, grepl("^020201|^020202|^020203|^020204|^020205|^020208|^02040|^02050|^020602", bnf_code))
drug_data_bnf$category <- "antihypertensive"

# first merge data.frame by BNF_Presentation_Code if the full code is available
drug_data_bnf <- merge(drug_data_bnf, antihyper_bnf_df[, c("BNF_Presentation_Code", "type", "drug")], all.x = T, by.x = "bnf_code", by.y = "BNF_Presentation_Code")

# identify the remaining drugs/types through string matching

types <- unique(antihyper_bnf_df$type)

for (t in types){
  antihyper_bnf_df_type <- antihyper_bnf_df[antihyper_bnf_df$type == t, ]
  code <- as.character(antihyper_bnf_df_type[ , "code"])[1]
  drugs <- unique(antihyper_bnf_df_type$drug)

  for (d in drugs){
    kws <- paste(unique(as.vector(antihyper_bnf_df[antihyper_bnf_df$drug == d, "root"])),collapse = '|')
    #print(kws)
    if (d == "sodium nitroprusside"){
      kws <- "nitroprusside"
    } else if (d == "diltiazem"){
      kws <- paste(c(kws, "adizem"), collapse = '|')
    } else if (d == "bendroflumethiazide"){
      kws <- paste(c(kws, "bendrofluazide"), collapse = '|')
    }

    drug_data_bnf <- drug_data_bnf %>% 
                      mutate(type = case_when(
                        !is.na(type) ~ type,
                        (is.na(type) & str_starts(bnf_code, code)) ~ t
                      ))

    drug_data_bnf <- drug_data_bnf %>% 
                      mutate(drug = case_when(
                        !is.na(drug) ~ drug,
                        (is.na(drug) & str_detect(drug_name, kws) & str_starts(bnf_code, code)) ~ d
                        ))    
  }
}

# do a string search per beta blocker drug # 02.04.0
bb_bnf_df = read.csv(snakemake@input[["bnf"]], sep = '\t')
# change type to betablocker
bb_bnf_df$type = "betablocker"
bb_bnf_df$keyword = tolower(bb_bnf_df$keyword)
drugs = as.vector(unique(bb_bnf_df$drug))

for (d in drugs){
  code = "02040"
  t <- as.character(bb_bnf_df[bb_bnf_df$drug == d, "type"])[1]
  kws = paste(as.vector(bb_bnf_df[bb_bnf_df$drug == d, "keyword"]),collapse = '|')

  drug_data_bnf <- drug_data_bnf %>% 
                      mutate(type = case_when(
                        !is.na(type) ~ type,
                        (is.na(type) & str_starts(bnf_code, code)) ~ t
                      ))

  drug_data_bnf <- drug_data_bnf %>% 
                      mutate(drug = case_when(
                        !is.na(drug) ~ drug,
                        (is.na(drug) & str_detect(drug_name, kws) & str_starts(bnf_code, code)) ~ d
                        )) 
}

# concatenate data and select columns of interest
cols = c("eid", "issue_date", "drug_name", "quantity", "drug", "type", "category")
drug_data <- rbind(drug_data_read2[, cols], drug_data_bnf[, cols])

# change category of beta-blocker
drug_data[(drug_data$type == "betablocker") & (!is.na(drug_data$type)), "category"] <- "betablocker"

# mark beta-blocker combination therapies 
combis = c("clopamide", "viskaldix", "hydchloroth", "sotazide", "tolerzide",
                "gppe", "bendroflumeth", "moducren", "prestim", "monozide", "tenben",
                "aspirin", "secadrex", "kalten", "spironol", "inderetic", "inderex", 
                "spiroprop", "co-betaloc", "lopresoretic", "corgaretic")
drug_data$bb_combi <- as.character(NA)

for (i in 1:length(combis)){
  c = combis[i]
  
  drug_data = drug_data %>% 
                    mutate(bb_combi = case_when(
                      ((str_detect(drug_name, c)) & (category == "betablocker")) ~ c,
                       !is.na(bb_combi) ~ bb_combi))
  i = i+1
}

#### manual corrections

## combinations with beta-blocker atenolol
drug_combis = c("beta-adalat", "co-tenidone", "kalten", "tenif", "tenoret", "atenolol+nifedipine", "atenolol+bendrofluazide",
                "atenolol+co-amilozide")

for (drug_combi in drug_combis){
  drug_data = drug_data %>% 
                mutate(drug = ifelse(str_detect(drug_name, drug_combi), "atenolol", as.character(drug)),
                        type = ifelse(str_detect(drug_name, drug_combi), "betablocker", as.character(type)),
                        category = ifelse(str_detect(drug_name, drug_combi), "betablocker", as.character(category)),
                        bb_combi = ifelse(str_detect(drug_name, drug_combi), drug_combi, as.character(bb_combi)))

}

## combinations with beta-blocker metoprolol and others
drug_combis = c("metoprolol tart+hydrochlorothiazide", "propranolol hcl+bendrofluazide", "sotalol hcl+hydrochlorothiazide")
drugs = c("metoprolol", "propranolol", "sotalol")

for (i in 1:length(drug_combis)){
  drug_combi = drug_combis[i]
  drug_ = drugs[i]
  drug_data = drug_data %>% 
                mutate(drug = ifelse(str_detect(drug_name, drug_combi), drug_, as.character(drug)),
                        type = ifelse(str_detect(drug_name, drug_combi), "betablocker", as.character(type)),
                        category = ifelse(str_detect(drug_name, drug_combi), "betablocker", as.character(category)),
                        bb_combi = ifelse(str_detect(drug_name, drug_combi), drug_combi, as.character(bb_combi)))
}

## take out urinary retention drugs
retention = paste(c("alfuzosin", "combodart", "tamsulosin", "flomax", "bazetham", "tabphyn", "stronazon", "tabphyn",
"petyme", "pamsva", "diffundox", "prosurin", "morvesin", "contiflo", "contiflo", "maxtron", "galebon", "pinexel", 
"tamurex", "cositam", "losinate", "faramsil", "flectone", "tamfrex", "solifenacin"), collapse = "|")
drug_data <- drug_data %>% filter(!(str_detect(drug_name, retention)))

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


