module load R/4.1.2
export R_LIBS=/scratch/90days/uqawal15/R_libraries/new_rlib_4.1.2
### preparing CNV list
R
library(dplyr)
library(readxl)
library(stringr)
options(scipen=999)
cnvs <- read_excel("/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/devCNVs/CNV_List.xlsx", sheet=1)
cnvs <- as.data.frame(cnvs)
cnvs$ID <- 1:nrow(cnvs)
cnvs$chrom <- sapply(strsplit(cnvs$Location, ":"), "[[", 1)
cnvs$chrom <- ifelse(cnvs$chrom == "X", 23, cnvs$chrom)
cnvs$Start <- sapply(strsplit(sapply(strsplit(cnvs$Location, ":"), "[[", 2), "-"), "[[", 1)
cnvs$End <- sapply(strsplit(sapply(strsplit(cnvs$Location, ":"), "[[", 2), "-"), "[[", 2)
dels <- cnvs %>% filter(Genotype=="Deletion")
dels <- dels %>% select(chrom, Start, End, ID)
write.table(dels, "devCNVs_dels.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
dups <- cnvs %>% filter(Genotype=="Duplication")
dups <- dups %>% select(chrom, Start, End, ID)
write.table(dups, "devCNVs_dups.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


### start burden analyses #############################################
Input="/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/overall_burden/rare"
WKDIR="/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/devCNVs"
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
#PNS -N devCNVs_dup
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=10:mem=100GB

p="/afm01/UQ/Q4399/software/plink-1.07-x86_64/plink"
WKDIR=/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/devCNVs

for i in {1..22}
do
region=$(awk 'NR == '${i}' {print $4}' ${WKDIR}/devCNVs_dups.txt)
echo ${region}
awk 'NR == '${i}' {print $1, $2, $3, $4}' ${WKDIR}/devCNVs_dups.txt > ${WKDIR}/region.txt
mkdir -p ${WKDIR}/DUP/${region}
${p} --noweb --cfile ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP --cnv-intersect ${WKDIR}/region.txt --cnv-overlap 0.5 --cnv-write --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}
${p} --noweb --cfile ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region} --cnv-make-map --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}

${p} --noweb --cfile ${WKDIR}/DUP/UKBB_CNVs_for_AN_hg38.DUP --cnv-intersect  ${WKDIR}/region.txt --cnv-region-overlap 0.5 --cnv-write --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}.reciprocal
${p} --noweb --cfile ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}.reciprocal --cnv-make-map --out ${WKDIR}/DUP/${region}/UKBB_CNVs_for_AN_hg38.DUP.${region}.reciprocal
done

#!/bin/bash
# 
#PBS -S /bin/bash
#PBS -A UQ-IMB
#PNS -N devCNVs_del
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=10:mem=100GB

p="/afm01/UQ/Q4399/software/plink-1.07-x86_64/plink"
WKDIR=/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/devCNVs


for i in {1..45}
do
region=$(awk 'NR == '${i}' {print $4}' ${WKDIR}/devCNVs_dels.txt)
echo ${region}
awk 'NR == '${i}' {print $1, $2, $3, $4}' ${WKDIR}/devCNVs_dels.txt > ${WKDIR}/region.txt
mkdir -p ${WKDIR}/DEL/${region}
${p} --noweb --cfile ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL --cnv-intersect ${WKDIR}/region.txt --cnv-overlap 0.5 --cnv-write --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}
${p} --noweb --cfile ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region} --cnv-make-map --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}

${p} --noweb --cfile ${WKDIR}/DEL/UKBB_CNVs_for_AN_hg38.DEL --cnv-intersect ${WKDIR}/region.txt --cnv-region-overlap 0.5 --cnv-write --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}.reciprocal
${p} --noweb --cfile ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}.reciprocal --cnv-make-map --out ${WKDIR}/DEL/${region}/UKBB_CNVs_for_AN_hg38.DEL.${region}.reciprocal
done

######################## changing the intersection rule for nrxn1
R
library(dplyr)
library(data.table)
df <- fread("/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/devCNVs/NRXN1_exons",fill=TRUE, sep=",") %>% as.data.frame()
exons <- data.frame()
for (i in 1:23){
  chr=2
  start=strsplit(df[1,"exonStarts"], split=",")[[1]][i]
  end=strsplit(df[1,"exonEnds"], split=",")[[1]][i]
  index=paste0(61, "_", i)
  vec <- c(chr, start, end, index)
  exons <- rbind(exons, vec)
}
write.table(exons, "NRXN1_exon_boundaries.bed", col.names=F, row.names=F, sep="\t", quote=F)


${p} --noweb --cfile ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61 --cnv-intersect ${WKDIR}/NRXN1_exon_boundaries.bed --cnv-write  --out ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.exons 
${p} --noweb --cfile ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.exons --cnv-make-map --out ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.exons 

${p} --noweb --cfile ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.reciprocal --cnv-intersect ${WKDIR}/NRXN1_exon_boundaries.bed --cnv-write  --out ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.reciprocal.exons 
${p} --noweb --cfile ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.reciprocal.exons --cnv-make-map --out ${WKDIR}/DEL/61/UKBB_CNVs_for_AN_hg38.DEL.61.reciprocal.exons 


## unweighted burden

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=devCNVs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /scratch/user/uqawal15/Anorexia/devCNVs.stdout
#SBATCH -e /scratch/user/uqawal15/Anorexia/devCNVs.stderr



module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

cd /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs

R --file=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs.R



###########################################################################
module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1
R
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

pos_dels <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs_dels.txt", header=FALSE)
pos_dups <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs_dups.txt", header=FALSE)
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
  if (index==61){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".exons.cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.exons.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dels) {
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
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
res$onesidedPvalue <- ifelse(res$Index%in%c(65,66) & res$OR<1, res$Pvalue/2, 
                             ifelse(!res$Index%in%c(65,66) & res$OR>1, res$Pvalue/2, 1-(res$Pvalue/2)))

write.table(res, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs_unweighted.txt", col.names=T, row.names=F, quote=F, sep="\t")



#### Repeat association analyses with BMI using linear regression
module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1


R
library(dplyr)
pheno <- read.table("/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat", header=TRUE)
fam <- read.table("/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38_BMI.fam")
AN <- read.table("/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.fam")
fam$V5 <- factor(fam$V5)
fam$Array <- pheno$Array[match(fam$V1, pheno$V2)]
fam$Sex <- pheno$Sex[match(fam$V1, pheno$V2)]
fam$Age <- pheno$Age[match(fam$V1, pheno$V2)]
fam$V6 <- ifelse(fam$V6==-9, NA, fam$V6)
fam$AN <- AN$V6[match(fam$V1, AN$V1)]

pos_dels <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs_dels.txt", header=FALSE)
pos_dups <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs_dups.txt", header=FALSE)
dels <- as.character(pos_dels$V4)
dups <- as.character(pos_dups$V4)

cnvs <- c(dels, dups)
res_cases <- vector("list", length(cnvs))
res_controls <- vector("list", length(cnvs))


## unweighted cases
for (i in 1:length(cnvs)){
  print(i)
  index <- cnvs[i]
  print(index)
  if (index==61){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".exons.cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.exons.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dels) {
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  }
  ## Cases
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
    names(vec) <- c("Index","Samples","N_CNVs","Avg_wCNV", "Avg_woCNV", "Beta", "CI", "Pvalue", "se")
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
    names(vec) <- c("Index","Samples","N_CNVs","Avg_wCNV", "Avg_woCNV", "Beta", "CI", "Pvalue", "se")
    res_cases[[i]] <- vec
  }
}

## unweighted controls
for (i in 1:length(cnvs)){
  print(i)
  index <- cnvs[i]
  print(index)
  if (index==61){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".exons.cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.exons.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dels) {
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DEL/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DEL.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
    cnv <- rbind(cnv1, cnv2) %>% distinct()
  } else if (index%in%dups){
    cnv1 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".cnv"), header=TRUE)
    cnv2 <- read.table(paste0("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/DUP/", cnvs[i], "/UKBB_CNVs_for_AN_hg38.DUP.", cnvs[i], ".reciprocal.cnv"), header=TRUE)
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
    names(vec) <- c("Index","Samples","N_CNVs","Avg_wCNV", "Avg_woCNV", "Beta", "CI", "Pvalue", "se")
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
    names(vec) <- c("Index","Samples","N_CNVs","Avg_wCNV", "Avg_woCNV", "Beta", "CI", "Pvalue", "se")
    res_controls[[i]] <- vec
    print(vec)
  }
}



res_controls <- do.call(rbind, res_controls)
res_controls <- res_controls %>% as.data.frame()
res_cases <- do.call(rbind, res_cases)
res_cases <- res_cases %>% as.data.frame()
write.table(res_controls, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/BMI/devCNVs_BMI_controls_UKB.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(res_cases, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/BMI/devCNVs_BMI_cases_UKB.txt", col.names=T, row.names=F, quote=F, sep="\t")



#### make cnv tracks for 3q29 duplications
p="/afm01/UQ/Q4399/software/plink-1.07-x86_64/plink"

${p} --noweb --cfile UKBB_CNVs_for_AN_hg38.DUP.10.unique --cnv-make-map --out UKBB_CNVs_for_AN_hg38.DUP.10.unique
${p} --noweb --cfile UKBB_CNVs_for_AN_hg38.DUP.10.unique --cnv-track --out UKBB_CNVs_for_AN_hg38.DUP.10.unique

