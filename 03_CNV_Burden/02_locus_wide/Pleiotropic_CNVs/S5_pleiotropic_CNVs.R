
### preparing CNV list
R
library(dplyr)
library(readxl)
library(stringr)
options(scipen=999)
collins <- read_excel("/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/drCNVs/1-s2.0-S0092867422007887-mmc4.xlsx", sheet=1)
collins <- as.data.frame(collins)
dels <- collins %>% filter(`CNV Type`=="DEL")
dels <- dels %>% select(Chrom, Start, End, `Segment ID`)
write.table(dels, "collins_dels.txt", sep="\t", col.names=F, row.names=F, quote=F)
dups <- collins %>% filter(`CNV Type`=="DUP")
dups <- dups %>% select(Chrom, Start, End, `Segment ID`)
write.table(dups, "collins_dups.txt", sep="\t", col.names=F, row.names=F, quote=F)



### start burden analyses #############################################
Input="/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/overall_burden/rare"
WKDIR="/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/drCNVs"
p="/afm01/UQ/Q4399/software/plink-1.07-x86_64/plink"
cd $WKDIR
mkdir -p DUP
mkdir -p DEL

${p} --noweb --cfile ${Input}/UKBB_CNVs_for_AN_hg38.rare --cnv-del --cnv-write  --out ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL
${p} --noweb --cfile ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL --cnv-make-map --out ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL

${p} --noweb --cfile ${Input}/UKBB_CNVs_for_AN_hg38.rare --cnv-dup --cnv-write  --out ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP
${p} --noweb --cfile ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP --cnv-make-map --out ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP



#!/bin/bash
# 
#PBS -S /bin/bash
#PBS -A UQ-IMB
#PNS -N collins_dup
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=10:mem=100GB

p="/afm01/UQ/Q4399/software/plink-1.07-x86_64/plink"
WKDIR=/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/drCNVs

for i in {1..101}
do
region=$(awk 'NR == '${i}' {print $4}' ${WKDIR}/collins_dups.txt)
echo ${region}
awk 'NR == '${i}' {print $1, $2, $3, $4}' ${WKDIR}/collins_dups.txt > ${WKDIR}/region.txt
mkdir -p ${WKDIR}/DUP/${region}
${p} --noweb --cfile ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP --cnv-kb 100 --cnv-intersect ${WKDIR}/region.txt --cnv-overlap 0.5 --cnv-write --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}
${p} --noweb --cfile ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region} --cnv-make-map --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}

${p} --noweb --cfile ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP --cnv-intersect ${WKDIR}/region.txt --cnv-region-overlap 0.25 --cnv-write --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}.reciprocal
${p} --noweb --cfile ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}.reciprocal --cnv-make-map --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}.reciprocal
done

#!/bin/bash
# 
#PBS -S /bin/bash
#PBS -A UQ-IMB
#PNS -N collins_del
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=10:mem=100GB

p="/afm01/UQ/Q4399/software/plink-1.07-x86_64/plink"
WKDIR=/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/drCNVs


for i in {1..77}
do
region=$(awk 'NR == '${i}' {print $4}' ${WKDIR}/collins_dels.txt)
echo ${region}
awk 'NR == '${i}' {print $1, $2, $3, $4}' ${WKDIR}/collins_dels.txt > ${WKDIR}/region.txt
mkdir -p ${WKDIR}/DEL/${region}
${p} --noweb --cfile ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL --cnv-kb 100 --cnv-intersect ${WKDIR}/region.txt --cnv-overlap 0.5 --cnv-write --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}
${p} --noweb --cfile ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region} --cnv-make-map --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}

${p} --noweb --cfile ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL --cnv-intersect ${WKDIR}/region.txt --cnv-region-overlap 0.25 --cnv-write --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}.reciprocal
${p} --noweb --cfile ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}.reciprocal --cnv-make-map --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}.reciprocal
done


## unweighted burden

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=drCNVs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /scratch/user/uqawal15/Anorexia/drCNVs.stdout
#SBATCH -e /scratch/user/uqawal15/Anorexia/drCNVs.stderr



module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

cd /QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs

R --file=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/drCNVs.R



###########################################################################

library(dplyr)
library(logistf)
library(data.table)
options(scipen=999)
pheno <- read.table("/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat", header=TRUE)
fam <- read.table("/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.fam")
fam$V5 <- factor(fam$V5)
fam$V6 <- factor(fam$V6)
fam$Array <- pheno$Array[match(fam$V1, pheno$V2)]
fam$Sex <- pheno$Sex[match(fam$V1, pheno$V2)]
fam$Age <- pheno$Age[match(fam$V1, pheno$V2)]

pos_dels <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/collins_dels.txt", header=FALSE)
pos_dups <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/collins_dups.txt", header=FALSE)
dels <- as.character(pos_dels$V4)
dups <- as.character(pos_dups$V4)

cnvs <- c(dels, dups)
res_firth <- vector("list", length(cnvs))
unadjusted_firth <- vector("list", length(cnvs))
logistic <- vector("list", length(cnvs))
unadjusted_logistic <- vector("list", length(cnvs))

for (i in 1:length(cnvs)){
  print(i)
  index <- cnvs[i]
  if (index%in%dels){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
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
  #  Firth Logistic Regression
  firth = logistf(V6 ~ Sex + Age + Array + cnv, data = fam)
  beta <- coef(firth)[5]
  se <- sqrt(diag(vcov(firth)))[5]
  OR <- signif(exp(coef(firth)[5]),3)
  dp2 <- nchar(OR)-2
  CI <- paste0("(", round(exp(confint(firth)[5,1]), dp2), "-", round(exp(confint(firth)[5,2]), dp2), ")")
  pval <- signif(firth$prob[5],dp2)
  vec <- c(index, "Firth", case_N, cont_N, signif(cont_se,2), cont_freqper10K, cont_CI, OR, CI, pval, beta, se)
  names(vec) <- c("Index", "Method", "Case_N", "Cont_N", "Cont_SE", "Cont_Freq_per10K", "Cont_CI", 
                  "OR", "CI", "Pvalue", "beta", "se")
  res_firth[[i]] <- vec
  # Unadjusted Firth Logistic Regression
  firth = logistf(V6 ~ cnv, data = fam)
  beta <- coef(firth)[2]
  se <- sqrt(diag(vcov(firth)))[2]
  OR <- signif(exp(coef(firth)[2]),3)
  dp2 <- nchar(OR)-2
  CI <- paste0("(", round(exp(confint(firth)[2,1]), dp2), "-", round(exp(confint(firth)[2,2]), dp2), ")")
  pval <- signif(firth$prob[2],dp2)
  vec <- c(index, "Unadjusted_Firth", case_N, cont_N, signif(cont_se,2), cont_freqper10K, cont_CI, OR, CI, pval, beta, se)
  names(vec) <- c("Index", "Method", "Case_N", "Cont_N", "Cont_SE", "Cont_Freq_per10K", "Cont_CI", 
                  "OR", "CI", "Pvalue", "beta", "se")
  unadjusted_firth[[i]] <- vec
  # Normal Logistic Regression
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
  names(vec) <- c("Index", "Method", "Case_N", "Cont_N", "Cont_SE", "Cont_Freq_per10K", "Cont_CI", 
                  "OR", "CI", "Pvalue", "beta", "se")
  logistic[[i]] <- vec
  # Unadjusted Logistic Regression
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
  names(vec) <- c("Index", "Method", "Case_N", "Cont_N", "Cont_SE", "Cont_Freq_per10K", "Cont_CI", 
                  "OR", "CI", "Pvalue", "beta", "se")
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

write.table(res, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/drCNVs_unweighted.txt", col.names=T, row.names=F, quote=F, sep="\t")



