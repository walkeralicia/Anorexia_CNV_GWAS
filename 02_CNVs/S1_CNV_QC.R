
#========================================================================================================
# This script filters acquired UKB CNV calls with extra CNV-level QC, selects for CNVs within UKB AN
# cases and controls, and formats CNV data into plink cfile format (.fam and .cnv). See S2_CNV_QC.sh to
# generate the corresponding .map cfile. 


# Read in required R libraries ==========================================================================
library(dplyr)

# Filter the UKB CNV summary statistics=================================================================
## Individual samples were excluded if they had
## <30 CNVs or >200 CNVs
## 1 CNV > 1Mb long
## waviness factor >0.03 or <−0.03
## call rate <96% or
## LRR SD >0.35
stats <- read.table("/QRISdata/Q2909/pheno/RunID_670814_CNV/Files_for_Retman/Summary_statistics.dat", header = TRUE, fill = TRUE)
filtered_stats <- stats %>%
  filter(
    NumCNV_per_person < 30,
    between(WavinessFactor, -0.03, 0.03),
    call_rate >= 96,
    LRR_SD <= 0.35,
    NumCNV_per_person < 200
  )

# Read in UKB AN case and control IDs=========================================================================
cases <- read.table("/QRISdata/Q4399/Anorexia/UKB/AN_IDs/case_ids.txt", header=FALSE)
controls <- read.table("/QRISdata/Q4399/Anorexia/UKB/AN_IDs/controls_ids.txt", header=FALSE)

# Filter for UKB AN cases and controls that have CNV calls=====================================================

## UKB ID file to map IDs between CNV calls and other UKB phenotypes
map <- read.table("/QRISdata/Q2909/pheno/RunID_670814_CNV/ukb12505bridge14421.txt")

## cases (1,260 people)
cases2 <- cases %>%
  mutate(cnvID = map$V2[match(V1, map$V1)],
         passed = as.integer(cnvID %in% filtered_stats$f.eid)) %>%
  filter(passed == 1) %>%
  select(V1, cnvID)
write.table(cases2, "/QRISdata/Q4399/Anorexia/UKB/sample_selection/case_ids.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## controls (385,930 people)
controls2 <- controls %>%
  mutate(cnvID = map$V2[match(V1, map$V1)],
         passed = as.integer(cnvID %in% filtered_stats$f.eid)) %>%
  filter(passed == 1) %>%
  select(V1, cnvID)
write.table(controls2, "/QRISdata/Q4399/Anorexia/UKB/sample_selection/controls_ids.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


# Extract UKB phenotypes (Age, sex, SNP array) for AN cases and controls======================================

# 21003: Age when attended assessment centre
# 31:    Sex 0 = female, 1 = male
# 22000: Genotype measurement batch

pheno <- read.csv("/QRISdata/Q2909/pheno/RAP/Cognitive_120122.csv")
pheno2 <- pheno %>%
  select(f.eid = eid, Age = X21003.0.0, Sex = X31.0.0, Array = X22000.0.0) %>%
  mutate(Array = ifelse(Array < 1, "BiLEVE", "Axiom"))

cases2 <- cases2 %>% mutate(AN = 2)
controls2 <- controls2 %>% mutate(AN = 1)
all <- bind_rows(cases2, controls2)

pheno3 <- merge(all, pheno2, by.x="V1", by.y="f.eid", all.x=TRUE)
write.table(pheno3, "/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat", col.names=T, row.names=F, quote=F, sep="\t")


# Filter individual UKB CNV calls to exclude the CNVs from samples that failed QC and remove CNVs with < 10 probes, confidence score < 10=====

# 1,930,027 CNV calls
cnvs <- read.table("/QRISdata/Q2909/pheno/RunID_670814_CNV/Files_for_Retman/All_CNVs_for_UKBB.dat", header = TRUE, fill = TRUE)
cnvs_filtered <- cnvs %>%
  filter(Probe >= 10, f.eid %in% pheno3$cnvID, Conf >= 10)
write.table(cnvs_filtered, "/QRISdata/Q4399/Anorexia/UKB/cnv_data/All_CNVs_for_UKBB_filtered.dat", col.names=T, row.names=F, quote=F, sep="\t")




# Create plink cnv cfiles (.fam and .cnv)=====================================================================================

pheno <- read.table("/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat", header=TRUE)
cnvs <- read.table("/QRISdata/Q4399/Anorexia/UKB/cnv_data/All_CNVs_for_UKBB_filtered.dat", header=TRUE,fill=TRUE)
# Merge datasets
cnvs_pheno <- left_join(pheno, cnvs, by = c("cnvID" = "f.eid"))

# .cnv file==================================================
#FID     Family ID
#IID     Individual ID
#CHR     Chromosome
#BP1     Start position (base-pair)
#BP2     End position (base-pair)
#TYPE    Type of variant, e.g. 0,1 or 3,4 copies
#SCORE   Confidence score associated with variant 
#SITES   Number of probes in the variant
cnv_file <- cnvs_pheno %>%
  select(FID = V1, IID = V1, CHR = chr, BP1 = start, BP2 = end, TYPE = Type, SCORE = Conf, SITES = Probe) %>%
  filter(!is.na(SCORE))
write.table(cnv_file, "/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN.cnv", col.names=T, row.names=F, quote=F, sep="\t")

# .fam file (387,190 people) ===============================
#Family ID
#Individual ID
#Paternal ID
#Maternal ID
#Sex (1=male; 2=female; other=unknown)
#Phenotype

fam_file <- cnvs_pheno[, c("V1", "V1","Age", "Sex", "AN")]
fam_file$Pat_ID <- 0
fam_file$Mat_ID <- 0
fam_file <- fam_file[, c("V1", "V1","Pat_ID", "Mat_ID", "Sex", "AN")]
fam_file_distinct <- distinct(fam_file)
fam_file_distinct$Sex <- ifelse(fam_file_distinct$Sex==0, 2, 
                                ifelse(fam_file_distinct$Sex==1, 1, NA))
fam_file_distinct2 <- fam_file_distinct[!is.na(fam_file_distinct$V1),]

fam_file <- cnvs_pheno %>%
  select(V1, Age, Sex, AN) %>%
  mutate(Pat_ID = 0, Mat_ID = 0, f.eid_2 = V1) %>%
  select(f.eid_2, V1, Pat_ID, Mat_ID, Sex, AN) %>%
  distinct() %>%
  mutate(Sex = ifelse(Sex == 0, 2, ifelse(Sex == 1, 1, NA))) %>%
  filter(!is.na(V1))

write.table(fam_file, "/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN.fam", col.names=F, row.names=F, quote=F, sep="\t")

