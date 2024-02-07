
#======This script calculates the burden of pleiotropic CNVs within Anorexia patients==================================

# load required R libraries
library(dplyr)
library(logistf)
library(data.table)
options(scipen=999)

# read in files ==================
pheno_path <- "/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat" ## path to phenotype file for covariates
fam <- read.table(pheno_path, header=T)
file_name <- "UKBB_CNVs_for_AN_hg38"

drCNV_path <- "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs"
pos_dels <- read.table(paste(drCNV_path, "collins_dels.txt", sep="/"),header=FALSE)
pos_dups <- read.table(paste(drCNV_path, "collins_dups.txt", sep="/"),header=FALSE)
dels <- as.character(pos_dels$V4)
dups <- as.character(pos_dups$V4)

cnvs <- c(dels, dups)
res_firth <- vector("list", length(cnvs))
unadjusted_firth <- vector("list", length(cnvs))
logistic <- vector("list", length(cnvs))
unadjusted_logistic <- vector("list", length(cnvs))
names <- c("Index", "Method", "Case_N", "Cont_N", "Cont_SE", "Cont_Freq_per10K", "Cont_CI", "OR", "CI", "Pvalue", "beta", "se")

# begin CNV burden analyses ========================
for (i in 1:length(cnvs)){
  print(i)
  index <- cnvs[i]
  if (index%in%dels) {
    cnv1 <- read.table(paste0(drCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0(drCNV_path, "/DEL/", cnvs[i], "/", file_name, ".DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0(drCNV_path, "/DUP/", cnvs[i], "/", file_name, ".DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0(drCNV_path, "/DUP/", cnvs[i], "/", file_name, ".DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  }
  if (is.na(cnv[1,1])) next
  fam$cnv <- as.factor(as.vector(ifelse(fam$V1%in%cnv$FID, "yes", "no")))
  fam_cnvs <- fam[fam$cnv=="yes",]
  case_N <- table(fam_cnvs$V6)[2]
  cont_N <- table(fam_cnvs$V6)[1]
  contTotal <- table(fam$V6)[1]
  cont_prop <- cont_N/contTotal
  cont_freqper10K <- signif(cont_prop*10000,3)
  dp1 <- nchar(cont_freqper10K)-2
  cont_se <- sqrt(cont_prop*(1-cont_prop)/contTotal)
  cont_LCI <- 10000*(cont_prop-(1.96*cont_se))
  cont_UCI <- 10000*(cont_prop+(1.96*cont_se))
  cont_CI <- paste0("(", round(cont_LCI, dp1), "-", round(cont_UCI, dp1), ")")
  
  #  Firth Logistic Regression =====================================
  firth = logistf(V6 ~ Sex + Age + Array + cnv, data = fam)
  beta <- coef(firth)[5]
  se <- sqrt(diag(vcov(firth)))[5]
  OR <- signif(exp(coef(firth)[5]),3)
  dp2 <- nchar(OR)-2
  CI <- paste0("(", round(exp(confint(firth)[5,1]), dp2), "-", round(exp(confint(firth)[5,2]), dp2), ")")
  pval <- signif(firth$prob[5],dp2)
  vec <- c(index, "Firth", case_N, cont_N, signif(cont_se,2), cont_freqper10K, cont_CI, OR, CI, pval, beta, se)
  names(vec) <- names
  res_firth[[i]] <- vec
  # Unadjusted Firth Logistic Regression ==========================================
  firth = logistf(V6 ~ cnv, data = fam)
  beta <- coef(firth)[2]
  se <- sqrt(diag(vcov(firth)))[2]
  OR <- signif(exp(coef(firth)[2]),3)
  dp2 <- nchar(OR)-2
  CI <- paste0("(", round(exp(confint(firth)[2,1]), dp2), "-", round(exp(confint(firth)[2,2]), dp2), ")")
  pval <- signif(firth$prob[2],dp2)
  vec <- c(index, "Unadjusted_Firth", case_N, cont_N, signif(cont_se,2), cont_freqper10K, cont_CI, OR, CI, pval, beta, se)
  names(vec) <- names
  unadjusted_firth[[i]] <- vec
 
  # Normal Logistic Regression ===================================================
  m <- glm(V6 ~ Sex + Age + Array + cnv, data = fam, family = "binomial")
  beta <- summary(m)$coefficients[5,1]
  se <- summary(m)$coefficients[5,2]
  OR <- signif(exp(summary(m)$coefficients[5,1]), 2)
  dp <- nchar(OR)-2
  LCI <- round(exp(summary(m)$coefficients[5,1]-1.96*exp(summary(m)$coefficients[5,2])), dp)
  HCI <- round(exp(summary(m)$coefficients[5,1]+1.96*exp(summary(m)$coefficients[5,2])), dp)
  CI <- paste0("(", LCI, "-", HCI, ")")
  pval <- signif(summary(m)$coefficients[5,4],2)
  vec <- c(index, "Logistic", case_N, cont_N, signif(cont_se,2), cont_freqper10K, cont_CI, OR, CI, pval, beta, se)
  names(vec) <- names
  logistic[[i]] <- vec
  # Unadjusted Logistic Regression ==============================================
  m <- glm(V6 ~ cnv, data = fam, family = "binomial")
  beta <- summary(m)$coefficients[2,1]
  se <- summary(m)$coefficients[2,2]
  OR <- signif(exp(summary(m)$coefficients[2,1]), 2)
  dp <- nchar(OR)-2
  LCI <- round(exp(summary(m)$coefficients[2,1]-1.96*exp(summary(m)$coefficients[2,2])), dp)
  HCI <- round(exp(summary(m)$coefficients[2,1]+1.96*exp(summary(m)$coefficients[2,2])), dp)
  CI <- paste0("(", LCI, "-", HCI, ")")
  pval <- signif(summary(m)$coefficients[2,4],2)
  vec <- c(index, "Unadjusted_Logistic", case_N, cont_N, signif(cont_se,2), cont_freqper10K, cont_CI, OR, CI, pval, beta, se)
  names(vec) <- names
  unadjusted_logistic[[i]] <- vec
  
}

res_firth2 <- do.call(rbind, res_firth)
res_firth2 <- res_firth2 %>% as.data.frame()
unadjusted_firth2 <- do.call(rbind, unadjusted_firth)
unadjusted_firth2 <- unadjusted_firth2 %>% as.data.frame()
logistic2 <- do.call(rbind, logistic)
logistic2 <- logistic2 %>% as.data.frame()
unadjusted_logistic2 <- do.call(rbind, unadjusted_logistic)
unadjusted_logistic2 <- unadjusted_logistic2 %>% as.data.frame()

res <- rbind(res_firth2, unadjusted_firth2, logistic2, unadjusted_logistic2)
res$OR <- as.numeric(as.character(res$OR))
res$Pvalue <- as.numeric(as.character(res$Pvalue))
res$onesidedPvalue <- ifelse(res$OR>1, res$Pvalue/2, 1-(res$Pvalue/2))

write.table(res, paste(drCNV_path, "UKBB_drCNVs_burden.txt", sep="/"), col.names=T, row.names=F, quote=F, sep="\t")



