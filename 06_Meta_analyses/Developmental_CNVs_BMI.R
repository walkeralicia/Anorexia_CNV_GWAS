
#============= This script meta-analyses the UKB and ANGI developmental CNV burden results with BMI as the outcome =============================

#== Read in required R libraries==
library(readxl)
library(stringr)
library(dplyr)
options(scipen=999)

#== Set paths ==
wkdir="/Anorexia/UKB/burden_analysis/devCNVs/BMI"
data_path<-"data" ## path to CNV_List.xlsx
ANGI_data<-"ANGI_data" ## path to folder with ANGI meta-analyses results

#=== AN controls ==============================================================================================================================

# read in data ==================
ukb <- read.table(paste(wkdir, "UKBB_devCNVs_BMI_controls.txt", sep="/"), header=T)
angi <- read.table(paste(ANGI_data, "ANGI_devCNVs_BMI_controls.txt", sep="/"), header=T)
stouffer <- merge(angi, ukb, by = c("Index"), all=TRUE)

# stouffer's method =================
stouffer$Zvalue.x <- (-1)*sign(stouffer$Beta.x)*qnorm(stouffer$Pvalue.x/2)
stouffer$Zvalue.y <- (-1)*sign(stouffer$Beta.y)*qnorm(stouffer$Pvalue.y/2)
stouffer$Stouffer_Zvalue <- (stouffer$Zvalue.x + stouffer$Zvalue.y)/sqrt(2)
stouffer$Stouffer_Pvalue <- 2*pnorm(-(abs(stouffer$Stouffer_Zvalue)))
stouffer$OR_meta <- signif((stouffer$Beta.x/(stouffer$se.x*stouffer$se.x) + stouffer$Beta.y/(stouffer$se.y*stouffer$se.y)) / (1/(stouffer$se.x*stouffer$se.x)+1/(stouffer$se.y*stouffer$se.y)),3)
stouffer$se_meta <- 1/sqrt(1/(stouffer$se.x*stouffer$se.x) +1/(stouffer$se.y*stouffer$se.y))
stouffer$dp <- nchar(stouffer$OR_meta)-2
stouffer$LCI_meta <- ifelse(stouffer$dp==-1, round((stouffer$OR_meta-1.96*stouffer$se_meta), 2), round((stouffer$OR_meta-1.96*stouffer$se_meta), stouffer$dp))
stouffer$HCI_meta <- ifelse(stouffer$dp==-1, round((stouffer$OR_meta+1.96*stouffer$se_meta), 2), round((stouffer$OR_meta+1.96*stouffer$se_meta), stouffer$dp))
stouffer$CI_meta <- paste0("(", stouffer$LCI_meta, "-", stouffer$HCI_meta, ")")


# annotate the information about the pleiotropic CNVs to the result file =====================

cnvs <- read_excel(paste0(data_path, "CNV_List.xlsx", sep="/"), sheet=1)
cnvs <- as.data.frame(cnvs)
cnvs$ID <- 1:nrow(cnvs)
comb <- merge(cnvs, stouffer, by.x = "ID", by.y="Index", all.x=TRUE)


res <- comb %>% 
  select(ID, CNV, Location, Genotype, Size,
         ANGI_Beta := Beta.x, ANGI_CI := CI.x, ANGI_Zvalue := Zvalue.x, ANGI_Pvalue := Pvalue.x,  
         ANGI_N_CNVs := N_CNVs.x, ANGI_Avg_wCNV := Avg_wCNV.x, ANGI_Avg_woCNV := Avg_woCNV.x,
         UKB_Beta := Beta.y, UKB_CI := CI.y, UKB_Zvalue := Zvalue.y, UKB_Pvalue := Pvalue.y, 
         UKB_N_CNVs := N_CNVs.y, UKB_Avg_wCNV := Avg_wCNV.y, UKB_Avg_woCNV := Avg_woCNV.y,
         Stouffer_Zvalue, Stouffer_Pvalue, Beta_OR := OR_meta, Meta_CI := CI_meta)

final <- res %>%
  mutate(across(c(ANGI_Zvalue, UKB_Zvalue, ANGI_Pvalue, UKB_Pvalue, Stouffer_Zvalue, Stouffer_Pvalue), ~signif(., 3))) %>%
  mutate(Bonf_Stouffer_Pvalue = Stouffer_Pvalue * 67) %>%
  arrange(Stouffer_Pvalue)


write.csv(final, paste(wkdir, "devCNVs_BMI_controls.csv", sep="/"),  row.names=F)
paper <- final %>% filter(Stouffer_Pvalue<0.05)
write.csv(paper, paste(wkdir, "devCNVs_BMI_controls_sig.csv", sep="/"),  row.names=F)

#========= AN cases===============================================================================================================================

# read in data ==================
ukb <- read.table(paste(wkdir, "UKBB_devCNVs_BMI_cases.txt", sep="/"), header=T)
angi <- read.table(paste(data_path, "ANGI_devCNVs_BMI_cases.txt", sep="/"), header=T)
stouffer <- merge(angi, ukb, by = c("Index"), all=TRUE)

# stouffer's method =================
stouffer$Zvalue.x <- (-1)*sign(stouffer$Beta.x)*qnorm(stouffer$Pvalue.x/2)
stouffer$Zvalue.y <- (-1)*sign(stouffer$Beta.y)*qnorm(stouffer$Pvalue.y/2)
stouffer$Stouffer_Zvalue <- (stouffer$Zvalue.x + stouffer$Zvalue.y)/sqrt(2)
stouffer$Stouffer_Pvalue <- 2*pnorm(-(abs(stouffer$Stouffer_Zvalue)))
stouffer$OR_meta <- signif((stouffer$Beta.x/(stouffer$se.x*stouffer$se.x) + stouffer$Beta.y/(stouffer$se.y*stouffer$se.y)) / (1/(stouffer$se.x*stouffer$se.x)+1/(stouffer$se.y*stouffer$se.y)),3)
stouffer$se_meta <- 1/sqrt(1/(stouffer$se.x*stouffer$se.x) +1/(stouffer$se.y*stouffer$se.y))
stouffer$dp <- nchar(stouffer$OR_meta)-2
stouffer$LCI_meta <- ifelse(stouffer$dp==-1, round((stouffer$OR_meta-1.96*stouffer$se_meta), 2), round((stouffer$OR_meta-1.96*stouffer$se_meta), stouffer$dp))
stouffer$HCI_meta <- ifelse(stouffer$dp==-1, round((stouffer$OR_meta+1.96*stouffer$se_meta), 2), round((stouffer$OR_meta+1.96*stouffer$se_meta), stouffer$dp))
stouffer$CI_meta <- paste0("(", stouffer$LCI_meta, "-", stouffer$HCI_meta, ")")


# annotate the information about the pleiotropic CNVs to the result file =====================

comb <- merge(cnvs, stouffer, by.x = "ID", by.y="Index", all.x=TRUE)

res <- comb %>% 
  select(ID, CNV, Location, Genotype, Size, 
         Beta.x, CI.x, Zvalue.x, Pvalue.x,  
         N_CNVs.x, Avg_wCNV.x, Avg_woCNV.x,
         Beta.y, CI.y, Zvalue.y, Pvalue.y, 
         N_CNVs.y, Avg_wCNV.y, Avg_woCNV.y,
         Stouffer_Zvalue, Stouffer_Pvalue, OR_meta, CI_meta)


res <- comb %>% 
  select(ID, CNV, Location, Genotype, Size,
         ANGI_Beta := Beta.x, ANGI_CI := CI.x, ANGI_Zvalue := Zvalue.x, ANGI_Pvalue := Pvalue.x,  
         ANGI_N_CNVs := N_CNVs.x, ANGI_Avg_wCNV := Avg_wCNV.x, ANGI_Avg_woCNV := Avg_woCNV.x,
         UKB_Beta := Beta.y, UKB_CI := CI.y, UKB_Zvalue := Zvalue.y, UKB_Pvalue := Pvalue.y, 
         UKB_N_CNVs := N_CNVs.y, UKB_Avg_wCNV := Avg_wCNV.y, UKB_Avg_woCNV := Avg_woCNV.y,
         Stouffer_Zvalue, Stouffer_Pvalue, Beta_OR := OR_meta, Meta_CI := CI_meta)

final <- res %>%
  mutate(across(c(ANGI_Zvalue, UKB_Zvalue, ANGI_Pvalue, UKB_Pvalue, Stouffer_Zvalue, Stouffer_Pvalue), ~signif(., 3))) %>%
  mutate(Bonf_Stouffer_Pvalue = Stouffer_Pvalue * 67) %>%
  arrange(Stouffer_Pvalue)


write.csv(final, paste(wkdir, "devCNVs_BMI_cases.csv", sep="/"),  row.names=F)
paper <- final %>% filter(Stouffer_Pvalue<0.05)
write.csv(paper, paste(wkdir, "devCNVs_BMI_cases_sig.csv", sep="/"),  row.names=F)

