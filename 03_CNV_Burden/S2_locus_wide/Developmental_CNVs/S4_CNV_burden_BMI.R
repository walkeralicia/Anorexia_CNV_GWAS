
#=============This script conducts CNV burden analyses with BMI as the outcome (within AN cases and controls)================

#== Load required R libraries ==

library(dplyr)
options(scipen=999)

#== Set up file paths ==

pheno_path <- "/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat" ## path to phenotype file for covariates
file_name <- "UKBB_CNVs_for_AN_hg38" ## cfile name
devCNV_path <- "/Anorexia/UKB/burden_analysis/devCNVs" ## path to calculate devCNV burden in


#== Read in files ==
fam <- read.table(pheno_path, header=T)

pos_dels <- read.table(paste(devCNV_path, "devCNVs_dels.txt",sep="/"), header=FALSE)
pos_dups <- read.table(paste(devCNV_path, "devCNVs_dups.txt",sep="/"), header=FALSE)
dels <- as.character(pos_dels$V4)
dups <- as.character(pos_dups$V4)

cnvs <- c(dels, dups)
res_cases <- vector("list", length(cnvs))
res_controls <- vector("list", length(cnvs))
names <- c("Index","Samples","N_CNVs","Avg_wCNV", "Avg_woCNV", "Beta", "CI", "Pvalue", "se")


## Cases
for (i in 1:length(cnvs)){
  print(i)
  index <- cnvs[i]
  print(index)
  if (index==61){
    cnv1 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".exons.cnv"), header=TRUE)
    cnv2 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".reciprocal.exons.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dels) {
    cnv1 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0(devCNV_path, "/DUP/", cnvs[i], "/", file_name, ".DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0(devCNV_path, "/DUP/", cnvs[i], "/", file_name, ".DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  }
  fam$cnv <- as.factor(as.character(as.vector(ifelse(fam$V1%in%cnv$FID, "yes", "no"))))
  fam_case <- fam %>% filter(AN==2)
  cnv_yes <- fam_case %>% filter(cnv=="yes")
  if (is.na(cnv_yes[1,1])) next
  avg_cnv_yes <- signif(mean(cnv_yes$V6, na.rm=TRUE),3)
  if (is.nan(avg_cnv_yes)){
    cnv_no <- fam_case %>% filter(cnv=="no")
    avg_cnv_no <- signif(mean(cnv_no$V6, na.rm=TRUE),3)
    n_cnvs <- length(cnv_yes$V1)
    vec <- c(index,"Cases", n_cnvs, avg_cnv_yes, avg_cnv_no, NA, NA, NA, NA)
    names(vec) <- names
    res_cases[[i]] <- vec
  } else {
    m <- lm(V6 ~ Sex + Age + Array + cnv, data = fam_case)
    m_summary <- summary(m)
    beta <- signif(m_summary$coefficients[5,1],2)
    dp <- nchar(beta)-2
    se <- m_summary$coefficients[5,2]
    p <- signif(m_summary$coefficients[5,4],2)
    LCI <- round(beta-1.96*se, dp)
    UCI <- round(beta+1.96*se, dp)
    CI <- paste0("(", LCI, "-", UCI, ")")
    cnv_yes <- fam_case %>% filter(cnv=="yes")
    avg_cnv_yes <- signif(mean(cnv_yes$V6, na.rm=TRUE),3)
    cnv_no <- fam_case %>% filter(cnv=="no")
    avg_cnv_no <- signif(mean(cnv_no$V6, na.rm=TRUE),3)
    n_cnvs <- length(cnv_yes$V1)
    vec <- c(index, "Cases", n_cnvs, avg_cnv_yes, avg_cnv_no, beta, CI, p, se)
    names(vec) <- names
    res_cases[[i]] <- vec
  }
}

## Controls
for (i in 1:length(cnvs)){
  print(i)
  index <- cnvs[i]
  print(index)
  if (index==61){
    cnv1 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".exons.cnv"), header=TRUE)
    cnv2 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".reciprocal.exons.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dels) {
    cnv1 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0(devCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0(devCNV_path, "/DUP/", cnvs[i], "/", file_name, ".DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0(devCNV_path, "/DUP/", cnvs[i], "/", file_name, ".DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  }
  fam$cnv <- as.factor(as.character(as.vector(ifelse(fam$V1%in%cnv$FID, "yes", "no"))))
  fam_cont <- fam %>% filter(AN==1)
  cnv_yes <- fam_cont %>% filter(cnv=="yes")
  if (is.na(cnv_yes[1,1])) next
  avg_cnv_yes <- signif(mean(cnv_yes$V6, na.rm=TRUE),3)
  if (is.nan(avg_cnv_yes)){
    cnv_no <- fam_cont %>% filter(cnv=="no")
    avg_cnv_no <- signif(mean(cnv_no$V6, na.rm=TRUE),3)
    n_cnvs <- length(cnv_yes$V1)
    vec <- c(index,"Controls", n_cnvs, avg_cnv_yes, avg_cnv_no, NA, NA, NA, NA)
    names(vec) <- names
    res_cases[[i]] <- vec
    print(vec)
  } else {
    m <- lm(V6 ~ Sex + Age + Array + cnv, data = fam_cont)
    m_summary <- summary(m)
    beta <- signif(m_summary$coefficients[5,1],2)
    dp <- nchar(beta)-2
    se <- m_summary$coefficients[5,2]
    p <- signif(m_summary$coefficients[5,4],2)
    if (dp==-1){
      LCI <- round(beta-1.96*se, 2)
      UCI <- round(beta+1.96*se, 2)
    } else {
      LCI <- round(beta-1.96*se, dp)
      UCI <- round(beta+1.96*se, dp)
    }
    CI <- paste0("(", LCI, "-", UCI, ")")
    cnv_yes <- fam_cont %>% filter(cnv=="yes")
    avg_cnv_yes <- signif(mean(cnv_yes$V6, na.rm=TRUE),3)
    cnv_no <- fam_cont %>% filter(cnv=="no")
    avg_cnv_no <- signif(mean(cnv_no$V6, na.rm=TRUE),3)
    n_cnvs <- length(cnv_yes$V1)
    vec <- c(index, "Controls", n_cnvs, avg_cnv_yes, avg_cnv_no, beta, CI, p, se)
    names(vec) <- names
    res_controls[[i]] <- vec
    print(vec)
  }
}


res_controls <- do.call(rbind, res_controls)
res_controls <- res_controls %>% as.data.frame()
res_cases <- do.call(rbind, res_cases)
res_cases <- res_cases %>% as.data.frame()
write.table(res_controls, paste(devCNV_path, "BMI/UKBB_devCNVs_BMI_controls.txt", sep="/"), col.names=T, row.names=F, quote=F, sep="\t")
write.table(res_cases, paste(devCNV_path, "BMI/UKBB_devCNVs_BMI_cases_UKB.txt", sep="/"), col.names=T, row.names=F, quote=F, sep="\t")




