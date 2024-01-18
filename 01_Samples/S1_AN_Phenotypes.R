

# This script extracts UKB individuals with the AN phenotype

# After running the script:
# Number of AN from each category
# 1. MHQ: 891 unique people
# 2. HESIN: ICD9_AN =1 record, ICD10_AN = 302 records, ICD10_AAN = 13 records (combined 100 unique people).
# 3. GP data: read3 = 247 unique people, read2 = 92 unique people. from read3 and read2, 339 unique people have AN (col=GP_anorexia_nervosa_case).
# 4. death due to AN: 5
# 5. only GP data that needs exclusion criteria (350 people) and cross validation (299 people) =  645 unique people with an overlap of 5. 

# Exclusions for category 5:
## 255 have have HES exclusion criteria and need to be excluded (col:anorexia_excluded)
## HES exclusion criteria is: ICD10_malignant_neo == 1 | other_exclusion_codes == 1 | malignant_neoplasms_ICD9 == 1 | neoplasms_ICD9 == 1 ~ 1
## 389 people leftover to possibly keep as cases (col:anorexia_not_excluded)
## A further 179 did not pass cross-validation with visiting eating disorder clinic (col = cross_val_excluded)
## 210 leftover after cross-validation 
# 6 people have diabetes

# total AN= 1398 (excluding the 6 with diabetes)

library(dplyr)
library(data.table)
library(tidyverse)

## mental health questionnaire
MHQ <- read.table(
  "/QRISdata/Q2909/pheno/10743_12505_UKBiobank_mentalHealth_260917.tab",
  sep="\t", header=TRUE)

df <- MHQ %>%
  mutate(AN_MHQ = 
           case_when(f.20544.0.1 == "16" |
             f.20544.0.2 == "16" |
             f.20544.0.3 == "16" |
             f.20544.0.4 == "16" |
             f.20544.0.5 == "16" |
             f.20544.0.6 == "16" |
             f.20544.0.7 == "16" |
             f.20544.0.8 == "16" |
             f.20544.0.9 == "16" |
             f.20544.0.10 == "16" |
             f.20544.0.11 == "16" |
             f.20544.0.12 == "16" |
             f.20544.0.13 == "16" |
             f.20544.0.14 == "16" |
             f.20544.0.15 == "16" |
             f.20544.0.16 == "16" ~ 1
           ))

df %>%
  count(AN_MHQ)

MHQ_dat <- df %>%
  select(f.eid,
         AN_MHQ)

colnames(MHQ_dat)
nrow(MHQ_dat)

#Filter out non-cases (to make dataset smaller)
MHQ_dat_cases <- MHQ_dat %>%
  filter(AN_MHQ == 1)

colnames(MHQ_dat_cases) <- c("eid", "AN_MHQ")


## HES inpatient diagnoses dataset 
# Within the hospital inpatient dataset, each inpatient episode for a participant is stored 
# as a single record, i.e. a row of data, or as multiple contiguous records. 
hesin_diag_diskframe <- read.table(
  "/QRISdata/Q2909/pheno/HESIN/hesin_diag_21122022.txt",
  sep="\t", header=TRUE)

hesin_diag <- hesin_diag_diskframe %>%
  select(eid,
         diag_icd9,
         diag_icd10,
  ) %>%
  collect()

### ICD9 - WHO International Classification of Diseases

ICD9 <- read.table(
  "/QRISdata/Q2909/pheno/RunID_43210_Aug2020/Diagnosis/ICD9_diagnosis.tab",
  sep="\t", header=TRUE)


# Anorexia nervosa
hesin_diag <- hesin_diag %>%
  mutate(ICD9_AN =
           case_when(diag_icd9 == "3071" ~ 1
           ))


hesin_diag %>%
  count(ICD9_AN) # 1


### ICD10 - WHO International Classification of Diseases

# Anorexia nervosa
hesin_diag <- hesin_diag %>%
  mutate(ICD10_AN =
           case_when(diag_icd10 == "F500" ~ 1
           ))

hesin_diag %>%
  count(ICD10_AN)

# Atypical anorexia nervosa
hesin_diag <- hesin_diag %>%
  mutate(ICD10_AAN =
           case_when(diag_icd10 == "F501" ~ 1
           ))
hesin_diag %>%
  count(ICD10_AAN)  # 13

#HES new dataframe of each disorder anorexia nervosa}
# Anorexia nervosa
hesin_AN <- hesin_diag %>%
  filter(ICD9_AN == 1 |
           ICD10_AN == 1 
  )
hesin_AN$ICD_AN <- 1
hesin_AN <- hesin_AN %>%
  select(eid, 
         ICD_AN)
# check no. of rows
nrow(hesin_AN)
# check no. of unique IDs 
length(unique(hesin_AN$eid))
# Keep only unique IDs
hesin_AN <- hesin_AN %>% 
  distinct(eid, .keep_all = TRUE) 

## HES new dataframe of each disorder atypical anorexia nervosa
# Atypical anorexia nervosa
hesin_AAN <- hesin_diag %>%
  filter(ICD10_AAN == 1)
hesin_AAN$ICD10_AAN <- 1 
hesin_AAN <- hesin_AAN %>%
  select(eid, 
         ICD10_AAN)
# check no. of rows
nrow(hesin_AAN)
# check no. of unique IDs 
length(unique(hesin_AAN$eid))
# Keep only unique IDs
hesin_AAN <- hesin_AAN %>% 
  distinct(eid, .keep_all = TRUE) 


##HES merge by ID
# Merge by ID
hesin_merged_ID <- list(
  hesin_AN,
  hesin_AAN) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
hesin_merged_ID %>%
  colnames()
hesin_merged_ID %>%
  nrow()
# Check same as nrow
length(unique(hesin_merged_ID$eid))

rm(hesin_AN)
rm(hesin_AAN)

# HES Exclusion criteria (for AN)
## Malignant neoplasms ICD9

### Convert empty spaces to NA
hesin_diag <- mutate_all(hesin_diag, list(~na_if(.,"")))
### Malignant neoplasms - ICD9
hesin_diag_exclusion <- hesin_diag %>%
  mutate(malignant_neoplasms_ICD9 =
           case_when(
             (diag_icd9 > 1400 &
                diag_icd9 < 2100) |
               
               (diag_icd9 > 140 &
                  diag_icd9 < 210) ~ 1,
             
             is.na(diag_icd9) ~ NA_real_,
             
             TRUE ~ 0
           )
  )

hesin_diag_exclusion %>%
  count(malignant_neoplasms_ICD9)

# Create dataset of malignant neoplasm cases
hesin_malignant_neo_ICD9_case <- hesin_diag_exclusion %>%
  filter(malignant_neoplasms_ICD9 == 1)


hesin_malignant_neo_ICD9_case <- hesin_malignant_neo_ICD9_case %>%
  select(eid, 
         malignant_neoplasms_ICD9)

# Check no. of unique IDs
length(unique(hesin_malignant_neo_ICD9_case$eid))

# Filter out non-unique IDs
hesin_malignant_neo_ICD9_case <- hesin_malignant_neo_ICD9_case %>% 
  distinct(eid,
           .keep_all = TRUE) 

#Create dataset of people without malignant neoplasms 

hesin_malignant_neo_ICD9_not_case <- hesin_diag_exclusion %>%
  filter(malignant_neoplasms_ICD9 == 0) 
hesin_malignant_neo_ICD9_not_case$no_malignant_neoplasms_ICD9 <- 1
hesin_malignant_neo_ICD9_not_case <- hesin_malignant_neo_ICD9_not_case %>%
  select(eid, 
         no_malignant_neoplasms_ICD9)
# check no. of rows
nrow(hesin_malignant_neo_ICD9_not_case)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_malignant_neo_ICD9_not_case$eid))
# Filter out non-unique IDs
hesin_malignant_neo_ICD9_not_case <- hesin_malignant_neo_ICD9_not_case %>% 
  distinct(eid, .keep_all = TRUE) 


#Merge the two datasets
hesin_malignant_neo_ICD_9_ID <- list(
  hesin_malignant_neo_ICD9_case,
  hesin_malignant_neo_ICD9_not_case
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
hesin_malignant_neo_ICD_9_ID %>%
  nrow()
length(unique(hesin_malignant_neo_ICD_9_ID$eid))
hesin_malignant_neo_ICD_9_ID %>%
  count(malignant_neoplasms_ICD9)
hesin_malignant_neo_ICD_9_ID %>%
  count(no_malignant_neoplasms_ICD9)
# Now have a dataset with no dup IDs, and one column for malignant neoplasms (1 or NA) 
#and the other for no malignant neoplasms (1 or NA). You can only be a true "no malignant neoplasms" 
#if you also do NOT have a 1 in the other column. (Each row is a data entry, not a participant)


## Neoplasms ICD9
### HES data exclusion criteria ICD9 neoplasms
## Neoplasms - ICD9
hesin_diag_exclusion <- hesin_diag_exclusion %>%
  mutate(neoplasms_ICD9 =
           case_when(
             (diag_icd9 > 2300 &
                diag_icd9 < 2400) |
               
               (diag_icd9 > 230 &
                  diag_icd9 < 240) ~ 1,
             
             is.na(diag_icd9) ~ NA_real_,
             
             TRUE ~ 0
           )
  )

hesin_diag_exclusion %>%
  count(neoplasms_ICD9) # 244

# Create dataset of cases
hesin_neo_ICD9_case <- hesin_diag_exclusion %>%
  filter(neoplasms_ICD9 == 1) 
hesin_neo_ICD9_case <- hesin_neo_ICD9_case %>%
  select(eid, 
         neoplasms_ICD9)
# check no. of rows
nrow(hesin_neo_ICD9_case)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_neo_ICD9_case$eid))
hesin_neo_ICD9_case <- hesin_neo_ICD9_case %>% 
  distinct(eid, .keep_all = TRUE) 

#HES data exclusion criteria NO ICD9 neoplasms
# Create dataset of people with NO neoplasms
hesin_neo_ICD9_not_case <- hesin_diag_exclusion %>%
  filter(neoplasms_ICD9 == 0) 
hesin_neo_ICD9_not_case$no_neoplasms_ICD9 <- 1
hesin_neo_ICD9_not_case <- hesin_neo_ICD9_not_case %>%
  select(eid, 
         no_neoplasms_ICD9)
# check no. of rows
nrow(hesin_neo_ICD9_not_case)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_neo_ICD9_not_case$eid))
hesin_neo_ICD9_not_case <- hesin_neo_ICD9_not_case %>% 
  distinct(eid, .keep_all = TRUE) 


# HES data exclusion criteria ICD9 neoplasms merge data}
# Merge the two datasets
hesin_neo_ID <- list(
  hesin_neo_ICD9_case,
  hesin_neo_ICD9_not_case
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
hesin_neo_ID %>%
  nrow()
length(unique(hesin_neo_ID$eid))
hesin_neo_ID %>%
  count(neoplasms_ICD9)
hesin_neo_ID %>%
  count(no_neoplasms_ICD9)


## Malignant neoplasms ICD10
## C00 -D09 & D37 -D49
### HES data exclusion criteria ICD10 malignant neoplasms}
M_NEO <- c("C0",
           "C1",
           "C2",
           "C3",
           "C4",
           "C5",
           "C6",
           "C7",
           "C8",
           "C9",
           "D0",
           "D37",
           "D38",
           "D39",
           "D40",
           "D41",
           "D42",
           "D43",
           "D44",
           "D45",
           "D46",
           "D47",
           "D48",
           "D49"
)

# MALIGNANT NEOPLASMS
# Identify all data entries of an ICD10 code of malignant neoplasm
hesin_diag_exclusion_ICD10 <- hesin_diag_exclusion %>%
  filter(stringr::str_detect(diag_icd10,
                   paste(M_NEO, collapse="|")))

# Create new column to indicate malginant neoplasm according to ICD-10
hesin_diag_exclusion_ICD10$ICD10_malignant_neo <- 1
hesin_diag_exclusion_ICD10 %>%
  count(ICD10_malignant_neo)
length(unique(hesin_diag_exclusion_ICD10$eid))

# Select only unique IDs
hesin_diag_exclusion_ICD10_ID <- hesin_diag_exclusion_ICD10 %>% 
  distinct(eid, .keep_all = TRUE) 


#Identify all data entries NOT of an ICD10 code of malignant neoplasm
###HES data exclusion criteria ICD10 NO malignant neoplasms}
# NO MALIGNANT NEOPLASMS
hesin_diag_NOT_exclusion_ICD10 <- hesin_diag_exclusion %>%
  filter(!str_detect(diag_icd10,
                     paste(M_NEO, collapse="|")))
# Test what happens to NAs
test_NA <- hesin_diag_NOT_exclusion_ICD10 %>%
  filter(is.na(diag_icd10))
test_NA %>% # No NAs
  nrow()
# Create new column to indicate NO malginant neoplasm according to ICD-10
hesin_diag_NOT_exclusion_ICD10$ICD10_NO_malignant_neo <- 1
hesin_diag_NOT_exclusion_ICD10 %>%
  count(ICD10_NO_malignant_neo)
# See how many of the data entries are from unique IDs
length(unique(hesin_diag_NOT_exclusion_ICD10$eid))
# Select only unique IDs
hesin_diag_NOT_exclusion_ICD10_ID <- hesin_diag_NOT_exclusion_ICD10 %>% 
  distinct(eid, .keep_all = TRUE) 


#Merge two datasets
##HES data exclusion criteria ICD10 malignant neoplasms merge datasets}
# Merge the two datasets
hesin_malignant_neo_ICD_10_ID <- list(
  hesin_diag_exclusion_ICD10_ID, # Participants who at some point have met exclusion criteria
  hesin_diag_NOT_exclusion_ICD10_ID # Participants who at some point have not met exclusion criteria (might overlap with above)
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Select only needed columns
hesin_malignant_neo_ICD_10_ID <- hesin_malignant_neo_ICD_10_ID %>%
  select(eid,
         ICD10_NO_malignant_neo,
         ICD10_malignant_neo)
# Check
hesin_malignant_neo_ICD_10_ID %>%
  nrow()
length(unique(hesin_malignant_neo_ICD_10_ID$eid)) # Number of unique IDs is same as number of rows = no dup IDs!
hesin_malignant_neo_ICD_10_ID %>%
  count(ICD10_NO_malignant_neo)
hesin_malignant_neo_ICD_10_ID %>%
  count(ICD10_malignant_neo)


## Further exclusion criteria ICD9 & ICD10 (i.e., beyond neoplasms)
# HES data exclusion criteria further exclusion criteria} 
# Cases after exclusion
hesin_diag_exclusion <- hesin_diag_exclusion %>%
  mutate(other_exclusion_codes =
           case_when(
             
             # ICD9 exclusion codes
             (diag_icd9 == "5645" |           # Functional diarrhea
                diag_icd9 == "2780" |            # Obesity
                diag_icd9 == "5569" |            # Idiopathic proctocolitis 
                diag_icd9 == "5303" |            # Stricture and stenosis of esophagus
                diag_icd9 == "5305" |            # Dyskinesia of esophagus
                diag_icd9 == "5308" |            # Other specified disorders of esophagus
                
                # ICD10 exclusion codes
                diag_icd10 == "K591"  |         # Functional diarrhea
                diag_icd10 == "E660"  |         # Obesity due to excess calories
                diag_icd10 == "K518"  |         # Other ulcerative colitis
                diag_icd10 == "K512"  |         # Ulcerative (chronic) proctitis
                diag_icd10 == "K513"  |         # Ulcerative (chronic) rectosigmoiditis
                diag_icd10 == "K514"  |         # Pseudopolyposis of colon
                diag_icd10 == "K515"  |         # Mucosal proctocolitis
                diag_icd10 == "K519"  |         # Ulcerative colitis, unspecified
                diag_icd10 == "K222"  |         # Oesophageal obstruction
                diag_icd10 == "K224"  |         # Dyskinesia of oesophagus
                diag_icd10 == "K228"  |         # Other specified diseases of oesophagus
                diag_icd10 == "K219"  |         # Gastro-oesophageal reflux disease without oesophagitis
                diag_icd10 == "J860"  |         # Pyothorax with fistula
                diag_icd10 == "K227"  |         # Barrett's oesophagus
                diag_icd10 == "K229"  |         # Disease of oesophagus, unspecified
                diag_icd10 == "F840"  |         # Childhood autism
                diag_icd10 == "F843"  |         # Other childhood disintegrative disorder
                diag_icd10 == "F845"  |         # Asperger's syndrome
                diag_icd10 == "F848"   |        # Other pervasive developmental disorders
                diag_icd10 == "F849")  ~ 1,     # Pervasive developmental disorder, unspecified
             
             diag_icd9 != "5645" ~ 0, # Every data entry that is NOT this one (and also would NOT be any of the others above due to the nature of case_when). This captures people who do have A entry in ICD9, but not one that meets exclusion criteria. (+++ HLD - to check this. I wrote this a while ago and now can't remember why I did it like this...)
             
             is.na(diag_icd9) |
               is.na(diag_icd10) ~ NA_real_, # Doesn't have any HES data (Entries with any non-exclusion criteria)
             
             TRUE ~ 0
             
           )
  )

# Check number
hesin_diag_exclusion %>% 
  count(other_exclusion_codes)
# Create dataset of cases
hesin_other_exclusion_codes_case <- hesin_diag_exclusion %>%
  filter(other_exclusion_codes == 1) 
hesin_other_exclusion_codes_case <- hesin_other_exclusion_codes_case %>%
  select(eid, 
         other_exclusion_codes)
# check no. of rows
nrow(hesin_other_exclusion_codes_case)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_other_exclusion_codes_case$eid))
hesin_other_exclusion_codes_case <- hesin_other_exclusion_codes_case %>% 
  distinct(eid, .keep_all = TRUE) 


#Create dataset of entries not meeting exclusion criteria 
##HES data exclusion criteria NO further exclusion criteria}
hesin_other_exclusion_codes_not_case <- hesin_diag_exclusion %>%
  filter(other_exclusion_codes == 0) 
hesin_other_exclusion_codes_not_case$no_other_exclusion_code <- 1
hesin_other_exclusion_codes_not_case <- hesin_other_exclusion_codes_not_case %>%
  select(eid, 
         no_other_exclusion_code)
# check no. of rows
nrow(hesin_other_exclusion_codes_not_case)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_other_exclusion_codes_not_case$eid))
hesin_other_exclusion_codes_not_case <- hesin_other_exclusion_codes_not_case %>% 
  distinct(eid, .keep_all = TRUE) 



#Merge the two datasets
## HES data exclusion criteria further exclusion criteria merge dataset}
hesin_other_codes_ID <- list(
  hesin_other_exclusion_codes_case,
  hesin_other_exclusion_codes_not_case
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
hesin_other_codes_ID %>%
  nrow()
length(unique(hesin_other_codes_ID$eid))
hesin_other_codes_ID %>%
  count(other_exclusion_codes)
hesin_other_codes_ID %>%
  count(no_other_exclusion_code)


### Merge all hesin exclusion data with 
###Merge exclusion datasets}
hesin_exclusion_dat <- list(
  hesin_malignant_neo_ICD_10_ID, # Malignant neoplasms ICD10
  hesin_other_codes_ID, # Other ICD9 or ICD10 codes
  hesin_neo_ID, # Neoplasms (ICD9)
  hesin_malignant_neo_ICD_9_ID, # Malignant neoplasms ICD9
  hesin_merged_ID # Dataset containing people who meet inclusion criteria via HES data
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
hesin_exclusion_dat %>%
  colnames()
hesin_exclusion_dat %>%
  nrow()
length(unique(hesin_exclusion_dat$eid))


##Identify people who meet at least one exclusion criteria
###identify people who meet at least one exclusion criteria}
hesin_exclusion_dat <- hesin_exclusion_dat %>%
  mutate(HES_exclusion_criteria =
           case_when(
             ICD10_malignant_neo == 1 |
               other_exclusion_codes == 1 |
               malignant_neoplasms_ICD9 == 1 |
               neoplasms_ICD9 == 1 ~ 1,
             
             TRUE ~ 0
           )
  )
hesin_exclusion_dat %>%
  count(HES_exclusion_criteria)

# Identify people with at least one HES data entry
hesin_exclusion_dat <- hesin_exclusion_dat %>%
  mutate(HES_data_entry =
           case_when(
             ICD10_NO_malignant_neo == 1 |
               no_other_exclusion_code == 1 |
               no_malignant_neoplasms_ICD9 == 1 |
               no_neoplasms_ICD9 == 1 ~ 1,
             
             TRUE ~ 0
           )
  )
hesin_exclusion_dat %>%
  count(HES_data_entry) # Most people have at least one HES data entry

#Remove dataframes not needed to save memory (only really need "hesin_exclusion_dat")
##HES remove dataframes to save memory}
rm(HES_data_cases)
rm(exclusion)
rm(hesin_diag)
rm(hesin_other_exclusion_codes_case)
rm(hesin_other_exclusion_codes_not_case)
rm(hesin_neo_ICD9_case)
rm(hesin_neo_ICD9_not_case)
rm(hesin_malignant_neo_ICD9_case)
rm(hesin_malignant_neo_ICD9_not_case)
rm(hesin_diag_exclusion_ICD10_ID)
rm(hesin_diag_NOT_exclusion_ICD10_ID)
rm(hesin_malignant_neo_ICD_10_ID)
rm(hesin_other_codes_ID)
rm(hesin_neo_ID)
rm(hesin_malignant_neo_ICD_9_ID)
gc()




## Death records: 
 
death_cause <- read.table(
  "/QRISdata/Q2909/pheno/DEATH/death_cause_221221.txt",
  sep="\t", header=TRUE)

# Death due to anorexia nervosa
death_data <- death_cause %>%
  mutate(death_due_to_AN =
           case_when(cause_icd10 == "F500" ~ 1
           ))

death_data %>%
  count(death_due_to_AN)
 
death_data_clean <- death_data %>%
  select(eid,
         death_due_to_AN)

death_data_cases <- death_data_clean %>%
  filter(death_due_to_AN == 1)

rm(death_data_clean)
rm(death_data)
rm(death_cause)

### Save data generated up until this point
dat_interim <- list(
  death_data_cases,
  hesin_exclusion_dat, # incl. info on exclusion cases
  MHQ_dat_cases
) %>% 
  reduce(full_join, 
         by = "eid"
  )
## check

nrow(dat_interim)
colnames(dat_interim)
length(unique(dat_interim$eid))
saveRDS(object = dat_interim, file = paste0("interim_death_MHQ_HES_UKB_PGCED3_case_data.rds"))

# Read in saved data 
dat_interim <- readRDS("interim_death_MHQ_HES_UKB_PGCED3_case_data.rds")



## GP clinical event records:
# Again, each row is a data entry rather than a participant.
gp_clinical_diskframe <- read.table(
  "/QRISdata/Q2909/pheno/GP/gp_clinical_211222.txt",
  sep="\t", header=TRUE, fill=TRUE)


gp_clinical <- gp_clinical_diskframe %>%
  select(eid,
         read_2,
         read_3,
         value1,
         value2,
         value3,
         event_dt # data of clinical event
         
  ) %>%
  collect()
head(gp_clinical)


## Anorexia nervosa
### GP read 2 code: Anorexia nervosa
### GP read code 2 anorexia nervosa}
# Eu500 = '[X]Anorexia nervosa'
gp_clinical <- gp_clinical %>%
  mutate(read2_Xanorexia_nervosa =
           case_when(read_2 == "Eu500" ~ 1
           )
  )
gp_clinical %>%
  count(read2_Xanorexia_nervosa)
# 1467. = 'H/OAnorexia nervosa'
gp_clinical <- gp_clinical %>%
  mutate(read2_HO_anorexia_nervosa =
           case_when(read_2 == "1467." ~ 1
           )
  )
gp_clinical %>%
  count(read2_HO_anorexia_nervosa)
# E271. = 'Anorexia nervosa'
gp_clinical <- gp_clinical %>%
  mutate(read2_anorexia_nervosa =
           case_when(read_2 == "E271." ~ 1
           )
  )

gp_clinical %>%
  count(read2_anorexia_nervosa)


#Create dataset
### GP read code 2 anorexia nervosa create dataset
gp_anorexia_nervosa_read2 <- gp_clinical %>%
  filter(read2_HO_anorexia_nervosa == 1 |
           read2_Xanorexia_nervosa == 1 |
           read2_anorexia_nervosa == 1) 
gp_anorexia_nervosa_read2$GP_anorexia_nervosa_read2 <- 1
gp_anorexia_nervosa_read2 <- gp_anorexia_nervosa_read2 %>%
  select(eid, 
         GP_anorexia_nervosa_read2)
# check no. of rows
nrow(gp_anorexia_nervosa_read2)
# check no. of unique IDs
length(unique(gp_anorexia_nervosa_read2$eid))
# Keep distinct IDs
gp_anorexia_nervosa_read2 <- gp_anorexia_nervosa_read2 %>% 
  distinct(eid, .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_anorexia_nervosa_read2 %>%
  nrow()
length(unique(gp_anorexia_nervosa_read2$eid))
# Check
gp_anorexia_nervosa_read2 %>%
  count(GP_anorexia_nervosa_read2)


### GP read 3 code: Anorexia nervosa 
### GP read code 3 anorexia nervosa}
# E271. = 'Anorexia nervosa'
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexia_nervosa =
           case_when(read_3 == "E271." ~ 1
           )
  )
gp_clinical %>%
  count(read3_anorexia_nervosa)
# 1467. = 'H/O: anorexia nervosa'
gp_clinical <- gp_clinical %>%
  mutate(read3_HO_anorexia_nervosa =
           case_when(read_3 == "1467." ~ 1
           )
  )
gp_clinical %>%
  count(read3_HO_anorexia_nervosa)
# .1467 = 'H/O: anorexia nervosa'
gp_clinical <- gp_clinical %>%
  mutate(read3_HO_anorexia_nervosa_2 =
           case_when(read_3 == ".1467" ~ 1
           )
  )
gp_clinical %>%
  count(read3_HO_anorexia_nervosa_2)
# .E49. = Anorexia nervosa and/or AN - Anorexia nervosa
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexia_nervosa_AN_anorexia_nervosa =
           case_when(read_3 == ".E49." ~ 1
           )
  )
gp_clinical %>%
  count(read3_anorexia_nervosa_AN_anorexia_nervosa)
# Eu500 = Anorexia nervosa and/or AN - Anorexia nervosa
gp_clinical <- gp_clinical %>%
  mutate(read3_Eu500_anorexia_nervosa_AN_anorexia_nervosa =
           case_when(read_3 == "Eu500" ~ 1
           )
  )
gp_clinical %>%
  count(read3_Eu500_anorexia_nervosa_AN_anorexia_nervosa)

##Create dataset
###GP read code 3 anorexia nervosa create dataset}
gp_anorexia_nervosa_read3 <- gp_clinical %>%
  filter(read3_anorexia_nervosa == 1  |
           read3_HO_anorexia_nervosa == 1 |
           read3_HO_anorexia_nervosa_2 == 1 |
           read3_anorexia_nervosa_AN_anorexia_nervosa == 1 |
           read3_Eu500_anorexia_nervosa_AN_anorexia_nervosa == 1
         
  ) 
gp_anorexia_nervosa_read3$GP_anorexia_nervosa_read3 <- 1
gp_anorexia_nervosa_read3 <- gp_anorexia_nervosa_read3 %>%
  select(eid, 
         GP_anorexia_nervosa_read3)
# check no. of rows
nrow(gp_anorexia_nervosa_read3)
# check no. of unique IDs 
length(unique(gp_anorexia_nervosa_read3$eid))
gp_anorexia_nervosa_read3 <- gp_anorexia_nervosa_read3 %>% 
  distinct(eid, .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_anorexia_nervosa_read3 %>%
  nrow()
length(unique(gp_anorexia_nervosa_read3$eid))
# Check
gp_anorexia_nervosa_read3 %>%
  count(GP_anorexia_nervosa_read3)


## Anorexia 
### GP read 2 code: Anorexia (these will need exclusion criteria)
## GP read code 2 anorexia need exclusion criteria}
# R030. = [D]Anorexia
gp_clinical <- gp_clinical %>%
  mutate(read2_DAnorexia =
           case_when(read_2 == "R030." ~ 1
           )
  )
gp_clinical %>%
  count(read2_DAnorexia)
# R030z = [D]Anorexia
gp_clinical <- gp_clinical %>%
  mutate(read2_AnorexiaNOS =
           case_when(read_2 == "R030z" ~ 1
           )
  )
gp_clinical %>%
  count(read2_AnorexiaNOS)


##Create dataset
###r GP read code 2 anorexia need create dataset
gp_anorexia_read2_need_excl_criteria <- gp_clinical %>%
  filter(read2_AnorexiaNOS == 1 |
           read2_DAnorexia == 1 
  ) 
gp_anorexia_read2_need_excl_criteria$GP_anorexia_read2_need_excl_criteria <- 1
gp_anorexia_read2_need_excl_criteria <- gp_anorexia_read2_need_excl_criteria %>%
  select(eid, 
         GP_anorexia_read2_need_excl_criteria)
# check no. of rows
nrow(gp_anorexia_read2_need_excl_criteria)
# check no. of unique IDs 
length(unique(gp_anorexia_read2_need_excl_criteria$eid))
gp_anorexia_read2_need_excl_criteria <- gp_anorexia_read2_need_excl_criteria %>% 
  distinct(eid,
           .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_anorexia_read2_need_excl_criteria %>%
  nrow()
length(unique(gp_anorexia_read2_need_excl_criteria$eid))
# Check
gp_anorexia_read2_need_excl_criteria %>%
  count(GP_anorexia_read2_need_excl_criteria) # 47


### GP read 3 code: Anorexia (these will need exclusion criteria)
## GP read code 3 anorexia need exclusion criteria}

# R030. = [D]Anorexia
gp_clinical <- gp_clinical %>%
  mutate(read3_Danorexia =
           case_when(read_3 == "R030." ~ 1
           )
  )
gp_clinical %>%
  count(read3_Danorexia)  # 98
# .R30. [D]Anorexia
gp_clinical <- gp_clinical %>%
  mutate(read3_Danorexia2 =
           case_when(read_3 == ".R30." ~ 1
           )
  )
gp_clinical %>%
  count(read3_Danorexia2)  # 0
# R030z = Anorexia NOS
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexiaNOS =
           case_when(read_3 == "R030z" ~ 1
           )
  )
gp_clinical %>%
  count(read3_anorexiaNOS)  # 37
#.R30Z [D]Anorexia NOS 
gp_clinical <- gp_clinical %>%
  mutate(read3_Danorexia_NOS =
           case_when(read_3 == ".R30Z" ~ 1
           )
  )
gp_clinical %>%
  count(read3_Danorexia_NOS)  # 0
# X76cG = Anorexia symptom
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexia_symptom =
           case_when(read_3 == "X76cG" ~ 1
           )
  )
gp_clinical %>%
  count(read3_anorexia_symptom)
# XE24f = Appetite loss - anorexia
gp_clinical <- gp_clinical %>%
  mutate(read3_appetite_loss_anorexia =
           case_when(read_3 == "XE24f" ~ 1
           )
  )
gp_clinical %>%
  count(read3_appetite_loss_anorexia)


##Create dataset
### GP read code 3 anorexia need create dataset
gp_anorexia_read3_need_excl_criteria <- gp_clinical %>%
  filter(read3_Danorexia == 1 |
           read3_Danorexia2 == 1 |
           read3_anorexiaNOS == 1 |
           read3_Danorexia_NOS == 1 |
           read3_anorexia_symptom == 1 |
           read3_appetite_loss_anorexia == 1
  ) 
gp_anorexia_read3_need_excl_criteria$GP_anorexia_read3_need_excl_criteria <- 1
gp_anorexia_read3_need_excl_criteria <- gp_anorexia_read3_need_excl_criteria %>%
  select(eid, 
         GP_anorexia_read3_need_excl_criteria)
# check no. of rows
nrow(gp_anorexia_read3_need_excl_criteria)
# check no. of unique IDs 
length(unique(gp_anorexia_read3_need_excl_criteria$eid))
gp_anorexia_read3_need_excl_criteria <- gp_anorexia_read3_need_excl_criteria %>% 
  distinct(eid,
           .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_anorexia_read3_need_excl_criteria %>%
  nrow()
length(unique(gp_anorexia_read3_need_excl_criteria$eid))
# Check
gp_anorexia_read3_need_excl_criteria %>%
  count(GP_anorexia_read3_need_excl_criteria)


### GP read 2 code: Anorexia (these will need exclusion criteria AND cross-validation)
#### Checking codes that are linked to multiple phenotypes
#For the below codes, we need more information to make a decision as to whether 
#to include them as 'anorexia nervosa' cases. This is because a single code is 
#linked to multiple phenotypes, some of which would be considered cases but others that would not be. 
#Specifically, some of these codes contain phenotypes that do not even mention anorexia, e.g., 'loss of appetite'

# 1612 = Appetite loss - anorexia OR Anorexia symptom 
gp_clinical %>%
  filter(read_2 == "1612.") %>%
  count(value1)
gp_clinical %>%
  filter(read_2 == "1612.") %>%
  count(value2)
gp_clinical %>%
  filter(read_2 == "1612.") %>%
  count(value3)

#Summary: Will need to cross-validate with attendance/referral to an eating disorder clinic.
#These cases need additional cross-validation as they are more likely to be anorexia in relation to cancer, 
#for example, because of the mention of 'loss of appetite' which is similar to the ICD10 code R360 which is 
#anorexia (not anorexia nervosa)

gp_clinical <- gp_clinical %>%
  mutate(read2_anorexia_appetite_loss_OR_anorexia_symptom =
           case_when(read_2 == "1612." ~ 1
           )
  )

gp_clinical %>%
  count(read2_anorexia_appetite_loss_OR_anorexia_symptom)


##Create dataset
### GP read code 2 anorexia need excl criteria and cross-validation create dataset
gp_anorexia_read2_need_excl_criteria_cross_val <- gp_clinical %>%
  filter(read2_anorexia_appetite_loss_OR_anorexia_symptom == 1 
  ) 
gp_anorexia_read2_need_excl_criteria_cross_val$GP_anorexia_read2_need_excl_criteria_cross_val <- 1
gp_anorexia_read2_need_excl_criteria_cross_val <- gp_anorexia_read2_need_excl_criteria_cross_val %>%
  select(eid, 
         GP_anorexia_read2_need_excl_criteria_cross_val)
# check no. of rows
nrow(gp_anorexia_read2_need_excl_criteria_cross_val)
# check no. of unique IDs 
length(unique(gp_anorexia_read2_need_excl_criteria_cross_val$eid))
gp_anorexia_read2_need_excl_criteria_cross_val <- gp_anorexia_read2_need_excl_criteria_cross_val %>% 
  distinct(eid,
           .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_anorexia_read2_need_excl_criteria_cross_val %>%
  nrow()
length(unique(gp_anorexia_read2_need_excl_criteria_cross_val$eid))
# Check
gp_anorexia_read2_need_excl_criteria_cross_val %>%
  count(GP_anorexia_read2_need_excl_criteria_cross_val)

### GP read 3 code: Anorexia (these will need exclusion criteria AND cross-validation with eating disorder clinic)
#### First, checking these codes that are linked to multiple phenotypes
#For the below codes, we need more information to make a decision as to whether to include them as 
#'anorexia nervosa' cases. This is because a single code is linked to multiple phenotypes, 
#'#some of which would be considered cases but others that would not be.

# XM07X = Lack of appetite OR Anorexia OR Off food OR No appetite OR Anorexic 
gp_clinical %>%
  filter(read_2 == "XM07X") %>%
  count(value1)
gp_clinical %>%
  filter(read_2 == "XM07X") %>%
  count(value2)
gp_clinical %>%
  filter(read_2 == "XM07X") %>%
  count(value3)
# 1612. = Appetite loss: [anorexia] or [anorexia symptom] OR Anorexia symptom OR Appetite loss - anorexia OR Loss of appetite OR Loss of appetite - symptom
gp_clinical %>%
  filter(read_2 == "1612.") %>%
  count(value1)
gp_clinical %>%
  filter(read_2 == "1612.") %>%
  count(value2)
gp_clinical %>%
  filter(read_2 == "1612.") %>%
  count(value3)
# .1612 = Anorexia symptom OR Appetite loss: [anorexia] or [anorexia symptom] OR Appetite loss - anorexia OR Loss of appetite OR Loss of appetite - symptom
gp_clinical %>%
  filter(read_2 == ".1612") %>%
  count(value1)
gp_clinical %>%
  filter(read_2 == ".1612") %>%
  count(value2)
gp_clinical %>%
  filter(read_2 == ".1612") %>%
  count(value3)


#Summary 18/02/22: No useful information. Will need to cross-validate with attendance/referral to an eating disorder clinic.

# XM07X = Lack of appetite OR Anorexia OR Off food OR No appetite OR Anorexic 
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexia_off_food =
           case_when(read_3 == "XM07X" ~ 1
           )
  )

gp_clinical %>%
  count(read3_anorexia_off_food)

# 1612. = Appetite loss: [anorexia] or [anorexia symptom] OR Anorexia symptom OR Appetite loss - anorexia OR Loss of appetite OR Loss of appetite - symptom
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexia_appetite_loss_OR_anorexia_symptom =
           case_when(read_3 == "1612." ~ 1
           )
  )


gp_clinical %>%
  count(read3_anorexia_appetite_loss_OR_anorexia_symptom)

# .1612 = Anorexia symptom OR Appetite loss: [anorexia] or [anorexia symptom] OR Appetite loss - anorexia OR Loss of appetite OR Loss of appetite - symptom
gp_clinical <- gp_clinical %>%
  mutate(read3_anorexia_appetite_loss_OR_anorexia_symptom2 =
           case_when(read_3 == ".1612" ~ 1
           )
  )


gp_clinical %>%
  count(read3_anorexia_appetite_loss_OR_anorexia_symptom2)



##Create dataset
###GP read code 3 anorexia need excl criteria and cross-validation create dataset
gp_anorexia_read3_need_excl_criteria_cross_val <- gp_clinical %>%
  filter(read3_anorexia_off_food == 1 |
           read3_anorexia_appetite_loss_OR_anorexia_symptom == 1 |
           read3_anorexia_appetite_loss_OR_anorexia_symptom2 == 1 
  ) 
gp_anorexia_read3_need_excl_criteria_cross_val$GP_anorexia_read3_need_excl_criteria_cross_val <- 1
gp_anorexia_read3_need_excl_criteria_cross_val <- gp_anorexia_read3_need_excl_criteria_cross_val %>%
  select(eid, 
         GP_anorexia_read3_need_excl_criteria_cross_val)
# check no. of rows
nrow(gp_anorexia_read3_need_excl_criteria_cross_val)
# check no. of unique IDs 
length(unique(gp_anorexia_read3_need_excl_criteria_cross_val$eid))
gp_anorexia_read3_need_excl_criteria_cross_val <- gp_anorexia_read3_need_excl_criteria_cross_val %>% 
  distinct(eid,
           .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_anorexia_read3_need_excl_criteria_cross_val %>%
  nrow()
length(unique(gp_anorexia_read3_need_excl_criteria_cross_val$eid))
# Check
gp_anorexia_read3_need_excl_criteria_cross_val %>%
  count(GP_anorexia_read3_need_excl_criteria_cross_val)


## Eating disorder clinic
### GP read 2 and 3 code: Seen in/referral to eating disorder clinics - used for cross-validation of cases we are unsure about
#i.e., People with a GP read 2 code of "1612." (Appetite loss - anorexia; Anorexia symptom;
#Loss of appetite - symptom), or read code 3 of "1612." (Appetite loss: [anorexia] or 
#[anorexia symptom]; Anorexia symptom; Appetite loss - anorexia; Loss of appetite; 
#Loss of appetite - symptom) or "XM07X" (Lack of appetite; Anorexia; Off food; No appetite; Anorexic) or 
#".1612" (Anorexia symptom; Appetite loss: [anorexia] or [anorexia symptom]; Appetite loss - anorexia; 
#Loss of appetite; Loss of appetite - symptom)

# 8HTN. = "Referral to eating disorders clinic"
gp_clinical <- gp_clinical %>%
  mutate(read2_referral_to_eating_disorder_clinic =
           case_when(read_2 == "8HTN." ~ 1
           )
  )
gp_clinical %>%
  count(read2_referral_to_eating_disorder_clinic)
# 9Nk9. Seen in eating disorder clinic
gp_clinical <- gp_clinical %>%
  mutate(read2_seen_in_eating_disorder_clinic =
           case_when(read_2 == "9Nk9." ~ 1
           )
  )
gp_clinical %>%
  count(read2_seen_in_eating_disorder_clinic)
# Xa2hW	Dietary advice for eating disorder
gp_clinical <- gp_clinical %>%
  mutate(read3_dietary_advice_for_eating_disorder =
           case_when(read_3 == "Xa2hW" ~ 1
           )
  )
gp_clinical %>%
  count(read3_dietary_advice_for_eating_disorder)
# XaEC2	Eating disorder counselling
gp_clinical <- gp_clinical %>%
  mutate(read3_eating_disorder_counselling =
           case_when(read_3 == "XaEC2" ~ 1
           )
  )
gp_clinical %>%
  count(read3_eating_disorder_counselling)
# XaPBE	Seen in eating disorder clinic
gp_clinical <- gp_clinical %>%
  mutate(read3_seen_in_eating_disorder_clinic =
           case_when(read_3 == "XaPBE" ~ 1
           )
  )
gp_clinical %>%
  count(read3_seen_in_eating_disorder_clinic)


##Create dataset
###GP read code 2 and 3 eating disorder clinic
gp_eating_disorder_clinic <- gp_clinical %>%
  filter(read3_seen_in_eating_disorder_clinic == 1 |
           read3_eating_disorder_counselling == 1 |
           read3_dietary_advice_for_eating_disorder == 1 |  
           read2_seen_in_eating_disorder_clinic == 1 |
           read2_referral_to_eating_disorder_clinic == 1
  ) 
gp_eating_disorder_clinic$GP_eating_disorder_clinic <- 1
gp_eating_disorder_clinic <- gp_eating_disorder_clinic %>%
  select(eid, 
         GP_eating_disorder_clinic)
# check no. of rows
nrow(gp_eating_disorder_clinic)
# check no. of unique IDs 
length(unique(gp_eating_disorder_clinic$eid))
gp_eating_disorder_clinic <- gp_eating_disorder_clinic %>% 
  distinct(eid,
           .keep_all = TRUE) 
# Check that same amount of unique IDs as there are rows
gp_eating_disorder_clinic %>%
  nrow()
length(unique(gp_eating_disorder_clinic$eid))
# Check
gp_eating_disorder_clinic %>%
  count(GP_eating_disorder_clinic)


### Merge all the GP AN data
###GP merge all AN data
# Merge the AN GP datasets
GP_dat_AN <- list(
  gp_anorexia_nervosa_read2,
  gp_anorexia_nervosa_read3,
  gp_anorexia_read2_need_excl_criteria,
  gp_anorexia_read3_need_excl_criteria,
  gp_anorexia_read2_need_excl_criteria_cross_val,
  gp_anorexia_read3_need_excl_criteria_cross_val,
  gp_eating_disorder_clinic
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
GP_dat_AN %>%
  colnames()
GP_dat_AN %>%
  nrow()


##Remove dataframes not needed
###interim remove dataframes to save memory
rm(gp_anorexia_nervosa_read2)                  
rm(gp_anorexia_nervosa_read3)
rm(gp_anorexia_read2_need_excl_criteria)
rm(gp_anorexia_read3_need_excl_criteria)
rm(gp_anorexia_read2_need_excl_criteria_cross_val)
rm(gp_anorexia_read3_need_excl_criteria_cross_val)
rm(gp_eating_disorder_clinic)
gc()


saveRDS(GP_dat_AN, "UKB_PGCED3_GP_AN_data.rds")
GP_dat_AN <- readRDS("UKB_PGCED3_GP_AN_data.rds")

# Merge ALL data: MHQ, Death, HES, and GP
dat <- list(
  GP_dat_AN,
  dat_interim
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
nrow(dat)
colnames(dat)
length(unique(dat$eid))
rm(dat_interim)

saveRDS(object = dat, paste0("interim2_death_MHQ_HES_GP_UKB_PGCED3_case_data.rds"))


# "Anorexia" and exclusion criteria
#CASES: "Anorexia" GP code and does NOT have HES exclusion code (incl. people without any HES data as cases), i.e., limited exclusions to people with known data on a counter indicated condition.


###anorexia exclusions prevalence of exlusionary conditions}
# Consider the prevalence of exclusionary conditions and see how likely it is that someone without HES data has an exclusion we don't know about
dat %>%
  count(HES_exclusion_criteria)
dat %>%
  count(HES_data_entry)
# Therefore, most people have at least one HES data entry and XX% of people with HES data meet exclusion criteria. 


#Create 'anorexia' variable (i.e., GP code of "anorexia") from both read codes that needs exclusion criteria, and one that needs exclusion criteria AND cross-validation
#anorexia exclusions create variable
# 'Anorexia' GP code needs excl criteria
dat <- dat %>%
  mutate(anorexia_excl_criteria =
           case_when(
             GP_anorexia_read2_need_excl_criteria == 1 |
               GP_anorexia_read3_need_excl_criteria == 1 ~ 1
           )
  )
# Check
dat %>%
  count(anorexia_excl_criteria)
# 'Anorexia' GP code needs excl criteria and cross-validation
dat <- dat %>%
  mutate(anorexia_excl_criteria_cross_val =
           case_when(
             GP_anorexia_read2_need_excl_criteria_cross_val == 1 |
               GP_anorexia_read3_need_excl_criteria_cross_val == 1 ~ 1
           )
  )
# Check
dat %>%
  count(anorexia_excl_criteria_cross_val)
# Check cross-over
dat %>%
  group_by(anorexia_excl_criteria) %>%
  count(anorexia_excl_criteria_cross_val)


##Do these people meet criteria for AN elsewhere? 
##anorexia exclusions meeting criteria elsewhere}
# Anorexia GP code AND does not meet criteria for anorexia nervosa elsewhere
dat <- dat %>%
  mutate(anorexia_GP_only_excl =
           case_when(
             anorexia_excl_criteria == 1 &
               (
                 is.na(AN_MHQ) &
                   is.na(death_due_to_AN) &
                   is.na(GP_anorexia_nervosa_read3) &
                   is.na(GP_anorexia_nervosa_read2) &
                   is.na(ICD_AN) &
                   is.na(ICD10_AAN)
               ) ~ 1
           )
  )
dat %>%
  count(anorexia_GP_only_excl)


# Anorexia GP code AND does not meet criteria for anorexia nervosa elsewhere
dat <- dat %>%
  mutate(anorexia_GP_only_excl_cross_val =
           case_when(
             anorexia_excl_criteria_cross_val == 1 &
               (
                 is.na(AN_MHQ) &
                   is.na(death_due_to_AN) &
                   is.na(GP_anorexia_nervosa_read3) &
                   is.na(GP_anorexia_nervosa_read2) &
                   is.na(ICD_AN) &
                   is.na(ICD10_AAN)
               ) ~ 1
           )
  )
dat %>%
  count(anorexia_GP_only_excl_cross_val)


###anorexia exclusions GP code only and meets exclusion criteria}
# EXCLUDE: Anorexia GP only and meets exclusion criteria
dat <- dat %>%
  mutate(anorexia_excluded =
           case_when(
             (anorexia_GP_only_excl == 1 |
                anorexia_GP_only_excl_cross_val == 1) &
               HES_exclusion_criteria == 1 ~ 1
           )
  )
dat %>%
  count(anorexia_excluded)
dat <- dat %>%
  mutate(anorexia_not_excluded =
           case_when(
             (anorexia_GP_only_excl == 1 |
                anorexia_GP_only_excl_cross_val == 1) &
               (HES_exclusion_criteria != 1 |
                  is.na(HES_exclusion_criteria)) ~ 1
           )
  )
dat %>%
  count(anorexia_not_excluded)



##Do these people have ANY HES data entry?
###anorexia exclusions HES data availability
dat %>%
  filter(anorexia_not_excluded == 1) %>%
  count(HES_data_entry) 


## Cross-validation step
###anorexia exclusions cross-validation
dat <- dat %>%
  mutate(anorexia_not_excluded_passed_cross_val =
           case_when(
             (anorexia_not_excluded == 1 & # Not yet excluded
                
                (anorexia_GP_only_excl_cross_val == 1 & # Needs cross-validation &...
                   is.na(anorexia_GP_only_excl)) & #...does not have ONLY the anorexia GP code (excl criteria only) and do not meet criteria for AN elsewhere
                
                GP_eating_disorder_clinic == 1 ~ 1), # passes cross-validation
             
             anorexia_not_excluded == 1 &
               anorexia_GP_only_excl == 1 ~ 1, # still need to carry forward the other cases with GP codes only needing exclusion critiera
             
           )
  )
# Check how many needing cross-validation pass cross-validation in the code above
dat %>%
  filter(anorexia_not_excluded == 1 & # Not yet excluded
           anorexia_GP_only_excl_cross_val == 1 & # Needs cross-validation &...
           is.na(anorexia_GP_only_excl) & #...does not have a GP code that only needs to pass exclusion criteria filter
           GP_eating_disorder_clinic == 1 )
# None passed cross-validation so have had to drop all of them
dat %>%
  count(anorexia_not_excluded_passed_cross_val)

dat_filt <- dat[!is.na(dat$anorexia_not_excluded) & is.na(dat$anorexia_not_excluded_passed_cross_val),]
dat$cross_val_excluded <- ifelse(dat$eid%in%dat_filt$eid, 1, 0)


## Additional exclusion criteria: 
##How many of these participants have type 1 diabetes? 
##additional exclusion criteria type 1 diabetes}
hesin_diag_diskframe <- read.table(
  "/QRISdata/Q2909/pheno/HESIN/hesin_diag_21122022.txt",
  sep="\t", header=TRUE)

hesin_diag <- hesin_diag_diskframe %>%
  select(eid,
         diag_icd9,
         diag_icd10,
  ) %>%
  collect()
nrow(hesin_diag)
colnames(hesin_diag)


##Diabetes check
##diabetes type 1 ICD9
# Type 1 diabetes ICD9
hesin_diag <- hesin_diag %>%
  mutate(type1_diabetes_ICD9 =
           case_when(diag_icd9 == "25001" |  # "Diabetes mellitus without mention of complication (juvenile type)"
                       diag_icd9 == "25011"   # "Diabetes with ketoacidosis (juvenile type)"
                     ~ 1
           ))
hesin_type1_diabetes_ICD9 <- hesin_diag %>%
  filter(type1_diabetes_ICD9 == 1) 
hesin_type1_diabetes_ICD9$type1_diabetes_ICD9 <- 1
hesin_type1_diabetes_ICD9 <- hesin_type1_diabetes_ICD9 %>%
  select(eid, 
         type1_diabetes_ICD9)
# check no. of rows
nrow(hesin_type1_diabetes_ICD9)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_type1_diabetes_ICD9$eid)) #17
hesin_type1_diabetes_ICD9 <- hesin_type1_diabetes_ICD9 %>% 
  distinct(eid, .keep_all = TRUE) 


##possible diabetes type 1 ICD9}
# Possible type 1 diabetes ICD9
hesin_diag <- hesin_diag %>%
  mutate(possible_type1_diabetes_ICD9 =
           case_when(
             diag_icd9 == "25009" | # "Diabetes mellitus without mention of compl. (adult/juvenile unspec.)"
               diag_icd9 == "25019" | # "Diabetes with ketoacidosis (adult/juvenile unspec.)"
               diag_icd9 == "25099" |  # "Diabetes with unspecified complications (unspecified onset)"
               diag_icd9 == "2503" | # "Diabetes with renal manifestations"
               diag_icd9 == "2504" | # "Diabetes with ophthalmic manifestations"
               diag_icd9 == "2505"  # "Diabetes with neurological manifestations"
             ~ 1
           ))
hesin_possible_type1_diabetes_ICD9 <- hesin_diag %>%
  filter(possible_type1_diabetes_ICD9 == 1) 
hesin_possible_type1_diabetes_ICD9$possible_type1_diabetes_ICD9 <- 1
hesin_possible_type1_diabetes_ICD9 <- hesin_possible_type1_diabetes_ICD9 %>%
  select(eid, 
         possible_type1_diabetes_ICD9)
# check no. of rows
nrow(hesin_possible_type1_diabetes_ICD9)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_possible_type1_diabetes_ICD9$eid))
hesin_possible_type1_diabetes_ICD9 <- hesin_possible_type1_diabetes_ICD9 %>% 
  distinct(eid, .keep_all = TRUE) 


## diabetes type 1 ICD10}
# Type 1 diabetes ICD10
hesin_diag <- hesin_diag %>%
  mutate(type1_diabetes_ICD10 =
           case_when( # Insulin-dependent diabetes mellitus...
             diag_icd10 == "E100" |  # With coma
               diag_icd10 == "E101" |   # With ketoacidosis
               diag_icd10 == "E102" |   # With renal complications
               diag_icd10 == "E103" |   # With ophthalmic complications
               diag_icd10 == "E104" |   # With neurological complications
               diag_icd10 == "E105" |   # With peripheral circulatory complications
               diag_icd10 == "E106" |   # With other specified complications
               diag_icd10 == "E107" |   # With multiple complications
               diag_icd10 == "E108" |   # With unspecified complications57
               diag_icd10 == "E109"   # Without complications
             ~ 1
           ))
hesin_type1_diabetes_ICD10 <- hesin_diag %>%
  filter(type1_diabetes_ICD10 == 1) 
hesin_type1_diabetes_ICD10$type1_diabetes_ICD10 <- 1
hesin_type1_diabetes_ICD10 <- hesin_type1_diabetes_ICD10 %>%
  select(eid, 
         type1_diabetes_ICD10)
# check no. of rows
nrow(hesin_type1_diabetes_ICD10)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_type1_diabetes_ICD10$eid)) # 4901
hesin_type1_diabetes_ICD10 <- hesin_type1_diabetes_ICD10 %>% 
  distinct(eid, .keep_all = TRUE) 
nrow(hesin_type1_diabetes_ICD10)


##possible diabetes type 1 ICD10}
# Possible type 1 diabetes ICD10
hesin_diag <- hesin_diag %>%
  mutate(possible_type1_diabetes_ICD10 =
           case_when(
             # Malnutrition-related diabetes mellitus...
             diag_icd10 == "E121" | # With ketoacidosis
               diag_icd10 == "E123" | # With ophthalmic complications
               diag_icd10 == "E125" | # With peripheral circulatory complications
               diag_icd10 == "E128"  | # With unspecified complications
               diag_icd10 == "E129"  | # Without complications
               
               # Other specified diabetes mellitus
               diag_icd10 == "E130" | # With ketoacidosis
               diag_icd10 == "E131" | # With ophthalmic complications
               diag_icd10 == "E132" | # With peripheral circulatory complications
               diag_icd10 == "E133"  | # With unspecified complications
               diag_icd10 == "E134"  | # Without complications
               diag_icd10 == "E135"  | # With peripheral circulatory complications
               diag_icd10 == "E136"  | # With other specified complications
               diag_icd10 == "E137"  | # With multiple complications
               diag_icd10 == "E138"  | # With unspecified complications 
               diag_icd10 == "E139" | # Without complications
               
               # Unspecified diabetes mellitus 
               diag_icd10 == "E140" | # With ketoacidosis
               diag_icd10 == "E141" | # With ophthalmic complications
               diag_icd10 == "E142" | # With peripheral circulatory complications
               diag_icd10 == "E143"  | # With unspecified complications
               diag_icd10 == "E144"  | # Without complications
               diag_icd10 == "E145"  | # With peripheral circulatory complications
               diag_icd10 == "E146"  | # With other specified complications
               diag_icd10 == "E147"  | # With multiple complications
               diag_icd10 == "E148"  | # With unspecified complications 
               diag_icd10 == "E149"  # Without complications
             ~ 1
           ))

hesin_possible_type1_diabetes_ICD10 <- hesin_diag %>%
  filter(possible_type1_diabetes_ICD10 == 1) 
hesin_possible_type1_diabetes_ICD10$possible_type1_diabetes_ICD10 <- 1
hesin_possible_type1_diabetes_ICD10 <- hesin_possible_type1_diabetes_ICD10 %>%
  select(eid, 
         possible_type1_diabetes_ICD10)
# check no. of rows
nrow(hesin_possible_type1_diabetes_ICD10)
# check no. of unique IDs - mismatch! need to get rid of dup IDs before merging
length(unique(hesin_possible_type1_diabetes_ICD10$eid))
hesin_possible_type1_diabetes_ICD10 <- hesin_possible_type1_diabetes_ICD10 %>% 
  distinct(eid, .keep_all = TRUE) 

###Merge all diabetes datasets by ID
##diabetes dat merge all datasets
diabetes_dat <- list(
  hesin_type1_diabetes_ICD9,
  hesin_possible_type1_diabetes_ICD9,
  hesin_type1_diabetes_ICD10,
  hesin_possible_type1_diabetes_ICD10
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
nrow(diabetes_dat) #9833
colnames(diabetes_dat)


##Drop unneeded datasets
##remove dataframes to save memory
rm(hesin_type1_diabetes_ICD9)
rm(hesin_possible_type1_diabetes_ICD9)
rm(hesin_type1_diabetes_ICD10)
rm(hesin_possible_type1_diabetes_ICD10)
gc()


### Merge diabetes data with data
##merge diabetes data with data
dat2 <- list(
  diabetes_dat,
  dat
) %>% 
  reduce(full_join, 
         by = "eid"
  )
# Check
nrow(dat2) #442440
colnames(dat2)


### Identify those of the 210 (all from anorexia GP code needing excl criteria only) with diabetes
##identify remaining anorexia cases with diabetes
dat2 <- dat2 %>%
  mutate(anorexia_not_excluded_diabetes =
           case_when(
             anorexia_not_excluded_passed_cross_val == 1 &
               (type1_diabetes_ICD9 == 1 |
                  type1_diabetes_ICD10 == 1) ~ 1
           )
  )
dat2 %>%
  count(anorexia_not_excluded_diabetes) # 2 definitely have type 1 diabetes
dat2 <- dat2 %>%
  mutate(anorexia_not_excluded_possible_diabetes =
           case_when(
             anorexia_not_excluded_passed_cross_val == 1 &
               (type1_diabetes_ICD9 == 1 |
                  type1_diabetes_ICD10 == 1 |
                  possible_type1_diabetes_ICD9 == 1 |
                  possible_type1_diabetes_ICD10 == 1) ~ 1
           )
  )
dat2 %>%
  count(anorexia_not_excluded_possible_diabetes) # an additional 6 possibly also have type 1 diabetes (3 of the people with definite type 1 diabetes also have possible diabetes - see code below)
# Check for overlap
dat2 %>%
  filter(anorexia_not_excluded_possible_diabetes == 1 &
           anorexia_not_excluded_diabetes == 1 ) %>%
  nrow()

# PGC numbers
## Anorexia nervosa cases
###PGC numbers anorexia nervosa overall cases}
# MHQ case
dat2 %>%
  count(AN_MHQ)

# GP case: Anorexia nervosa
dat2 <- dat2 %>%
  mutate(GP_anorexia_nervosa_case =
           case_when(   
             # GP data
             GP_anorexia_nervosa_read3 == 1 |
               GP_anorexia_nervosa_read2 == 1 ~ 1,
             
             TRUE ~ 0 # Not AN case
           )
  )
# Check
dat2 %>%
  count(GP_anorexia_nervosa_case)


# GP case: Anorexia but does not meet excl criteria (incl. diabetes)
dat2 <- dat2 %>%
  mutate(GP_anorexia_case =
           case_when(   
             # GP data
             (anorexia_not_excluded_passed_cross_val == 1 & # GP code of anorexia, does not meet AN criteria elsewhere, and does not meet original excl. criteria...
                (is.na(anorexia_not_excluded_diabetes) &
                   is.na(anorexia_not_excluded_possible_diabetes))) ~ 1, # ..and does not meet criteria for diabetes
             
             TRUE ~ 0 # Not AN case
           )
  )
# Check
dat2 %>%
  count(GP_anorexia_case) #people left over (additionally excluded 6 people of the 260 due to diabetes and thus possible diabulimia)


# HES case
dat2 <- dat2 %>%
  mutate(HES_AN_case =
           case_when(   
             # HES data
             ICD_AN == 1 | # ICD-10 code of F500 (Anorexia nervosa) or ICD-9 code of 3071 (Anorexia nervosa) 
               ICD10_AAN == 1 ~ 1, # or ICD-10 code of F501 (atypical Anorexia nervosa)
             
             TRUE ~ 0 # Not AN case
           )
  )
# Check
dat2 %>%
  count(HES_AN_case) # 100

# Overall AN cases
dat2 <- dat2 %>%
  mutate(AN_case =
           case_when( # MHQ data
             AN_MHQ == 1 | 
               
               # GP data
               GP_anorexia_nervosa_case == 1 |
               GP_anorexia_case == 1 |
               
               # Death data
               death_due_to_AN == 1 |
               
               # HES data
               HES_AN_case  ~ 1,
             
             TRUE ~ 0 # Not AN case
           )
  )


# Check
dat2 %>%
  count(AN_case) # 1404

#Explore where AN cases were identified
##PGC numbers anorexia nervosa origin case identification
table <- dat2 %>%
  filter(AN_case == 1) %>%
  group_by(GP_anorexia_nervosa_case,
           GP_anorexia_case,
           death_due_to_AN,
           AN_MHQ,
           HES_AN_case) %>%
  count()
table

# Double check how many are GP only (for info on additional cases obtained via medical records)
GP_cases <- dat2 %>% 
  filter(
    (
      (GP_anorexia_nervosa_case == 1 |
         GP_anorexia_case == 1) &
        
        ((death_due_to_AN == 0 |
            is.na(death_due_to_AN)) &
           (AN_MHQ == 0 |
              is.na(AN_MHQ)) &
           (HES_AN_case ==0 |
              is.na(HES_AN_case)))
    ) 
  )
nrow(GP_cases) # 436


# Save final dataset
##First all variables
##save data all variables
all_variables_case_dat <- dat2 %>%
  select(eid,
         
         # Cleaned analysis variables
         AN_case,
         
         # BE variables
         
         # AN variables
         GP_anorexia_nervosa_case,
         GP_anorexia_case,
         death_due_to_AN,
         AN_MHQ,
         HES_AN_case)
saveRDS(object = all_variables_case_dat,
        file = paste0("UKB_PGCED3_all_case_data_070322.rds"))


#Second, only cleaned variables necessary for analysis
##save data cleaned analysis variables only
only_case_dat <- dat2 %>%
  select(eid,
         AN_case) # 442440
saveRDS(object = only_case_dat, file = paste0("UKB_PGCED3_only_case_data_070322.rds"))



# create table of exclusion criteria:
exclusion_dat <- dat2 %>%
  select(eid,
         
         # Cleaned analysis variables
         AN_case,
         
         # HES exclusion criteria
         anorexia_excluded,
         # failed cross-validation
         cross_val_excluded,
         # diabetes but not excluded
         anorexia_not_excluded_diabetes,
         # possible diabeted but not excluded
         anorexia_not_excluded_possible_diabetes)
         #GP read 2 and read 3 exclusion criteria)
saveRDS(object = exclusion_dat,
        file = paste0("UKB_PGCED3_exclusion_dat.rds"))

exclusion_dat2 <- exclusion_dat[exclusion_dat$AN_case==1,]
dat_control_excluded <- exclusion_dat2[!is.na(exclusion_dat2$anorexia_excluded) |
                                         exclusion_dat2$cross_val_excluded==1 |
                                         !is.na(exclusion_dat2$anorexia_not_excluded_diabetes) |
                                         !is.na(exclusion_dat2$anorexia_not_excluded_possible_diabetes), ]

## covariates:
#Age (i.e., when assessed);
#Sex (biological sex at birth);
#Current BMI (at time of assessment)
#You need a file with required fields, one per line, no header.

#COVARIATES read in data
f <- read.table(
  "/QRISdata/Q2909/pheno/RunID_43210_Aug2020/Quant_heritable/Quantitative_pheno.tab",
  sep="\t", header=TRUE)


##COVARIATES add field code to file
strings <- c("f.eid",
             "f.31.0.0", # sex
             "f.21022.0.0", # age at recruitment
             "f.34.0.0", # Birth year
             "f.21001.0.0" # BMI (0 = Initial assessment visit)
) 

f2 <- f[, strings]

colnames(f2) <- c("eid", "sex", "age_at_recruitment", "birth_year", "bmi_at_initial_recruitment")
saveRDS(object = f2, file = paste0("small_phenotype_subset.rds"))



## PCs:
covar_PCs <- read.table(
  "/QRISdata/Q2909/v2Samples/ukbEUR_PC40_HM3excludeRegions_v2.txt",
  sep="\t", header=TRUE)

covar_PCs_select_columns <- covar_PCs %>%
  select(IID, PC1, PC2, PC3, PC4, PC5, PC6)

saveRDS(object = covar_PCs_select_columns, file = paste0("PC_covariates.rds"))


# Extract AN phenotype

## AN Cases that do not meet the inclusion criteria (440 people)
exclusion_dat <- readRDS("/QRISdata/Q4399/Anorexia/UKB/AN_IDs/UKB_PGCED3_exclusion_dat.rds")
controls <- exclusion_dat[exclusion_dat$AN_case==0,]
cases_excluded <- controls[!is.na(controls$anorexia_excluded) |
                             controls$cross_val_excluded==1 |
                             !is.na(controls$anorexia_not_excluded_diabetes) |
                             !is.na(controls$anorexia_not_excluded_possible_diabetes), ]
cases_excluded_ids <- cases_excluded$eid
write.table(cases_excluded_ids, "/QRISdata/Q4399/Anorexia/UKB/AN_IDs/cases_excluded_ids.txt", col.names=F, row.names=F, quote=F, sep="\t")

## AN cases (white and non-withdrawn) (1,327 people)
withdrawn <- read.table("/QRISdata/Q2909/v2Samples/ukb_sampleExclusions_v2_and_withdrawn_w12505_20220222.list")
eur <- read.table("/QRISdata/Q2909/new_ancestry_calls_2019/projection/EUR.id")

AN <- readRDS("/QRISdata/Q4399/Anorexia/UKB/AN_IDs/UKB_PGCED3_only_case_data_070322.rds")
AN_case <- AN[AN$AN_case==1,] # 1404
AN_case_check <- AN_case[!AN_case$eid%in%cases_excluded_ids,]
cases <- AN_case_check
cases$eur <- ifelse(cases$eid%in%eur$V1, 1, 0)
cases$withdrawn <- ifelse(cases$eid%in%withdrawn$V1, 1, 0)
cases2 <- cases[cases$eur==1 & cases$withdrawn == 0, ]
write.table(cases2, "/QRISdata/Q4399/Anorexia/UKB/AN_IDs/case_ids.txt", col.names=F, row.names=F, quote=F, sep="\t")


## AN controls (white and non-withdrawn) (406,256 people)
controls_keep <- controls[!controls$eid%in%cases_excluded_ids,]
controls <- controls_keep
controls$eur <- ifelse(controls$eid%in%eur$V1, 1, 0)
controls$withdrawn <- ifelse(controls$eid%in%withdrawn$V1, 1, 0)
controls2 <- controls[controls$eur==1 & controls$withdrawn == 0, ]
write.table(controls2, "/QRISdata/Q4399/Anorexia/UKB/AN_IDs/controls_ids.txt", col.names=F, row.names=F, quote=F, sep="\t")

