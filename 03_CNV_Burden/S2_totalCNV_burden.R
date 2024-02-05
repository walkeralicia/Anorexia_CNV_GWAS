
#=========== This script annotates rare CNVs and evaluates total CNV burden within the UKB============================================

# load required R packages
library(dplyr)
library(tidyverse)
library(data.table)
library(skimr)
library(rstatix)
library(stringr)
library(ggplot2)

#==== First extract collins et al. dosage sensitive genes ===========================

path_drCNVS="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs" ## path to disease-risk dosage-sensitive CNVs from Collins et al.


collins <- read.table(paste(path_drCNVs,"collins2022-TableS7-pHaplo-pTriplo.tsv", sep="/"), header=TRUE)
genes <- read.table(paste0(path_drCNVs, "geneMatrix.tsv"), header=TRUE, fill=TRUE)
genes2 <- genes[, c(1,8,9,10)]
both <- merge(genes2, collins, by="ensgid", all.y=TRUE)
haplo <- both %>% filter(pHaplo >= 0.86)
triplo <- both %>% filter(pTriplo >= 0.96)
haplo2 <- haplo %>% select(hg38h0, h1, h2, ensgid)
haplo2$hg38h0 <- gsub("chr", "", haplo2$hg38h0)
triplo2 <- triplo %>% select(hg38h0, h1, h2, ensgid)
triplo2$hg38h0 <- gsub("chr", "", triplo2$hg38h0)
write.table(triplo2, paste(path_drCNVs,"triplosensitive_list.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
write.table(haplo2, paste(path_drCNVs,"haploinsufficient_list.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")


#============ intersect rare UKB CNVs with disease-risk (dr) Collins dosage-sensitive genes ==================================

system("p=/QRISdata/Q4399/software/plink-1.07-x86_64/plink")
system("wkdir=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden")
system("drCNVS=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs")

system("${p} --noweb --cfile ${wkdir}/UKBB_CNVs_for_AN_hg38.DUP --cnv-intersect ${drCNVs}/triplosensitive_list.txt --cnv-write  --out ${wkdir}/UKBB_CNVs_for_AN_hg38.DUP.sensitive")
system("${p} --noweb --cfile ${wkdir}/UKBB_CNVs_for_AN_hg38.DUP.sensitive --cnv-make-map --out ${wkdir}/UKBB_CNVs_for_AN_hg38.DUP.sensitive")

system("${p} --noweb --cfile ${wkdir}/UKBB_CNVs_for_AN_hg38.DEL --cnv-intersect ${drCNVs}/haploinsufficient_list.txt --cnv-write  --out ${wkdir}/UKBB_CNVs_for_AN_hg38.DEL.sensitive")
system("${p} --noweb --cfile ${wkdir}/UKBB_CNVs_for_AN_hg38.DEL.sensitive --cnv-make-map --out ${wkdir}/UKBB_CNVs_for_AN_hg38.DEL.sensitive")

#================ Annotating rare CNVs with mammalian constraint (1K bins from Pat Sullivan)======================

wkdir="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden" ## path to conduct total CNV burden analyses in.

a <- fread(paste(wkdir, "UKBB_CNVs_for_AN_hg38.ALL.cnv", sep="/")) %>% 
  mutate(chr = ifelse(CHR==23, "chrX", paste0("chr", CHR)),
         tmp = row_number()) %>% 
  arrange(chr, BP1, BP2) %>% 
  select(chr, BP1, BP2, tmp)

f1 <- a %>% 
  mutate(chr = chr, 
         start0 = as.integer(1000*floor  (BP1/1000)),
         end    = as.integer(1000*ceiling(BP2/1000))) %>% 
  select(chr, start0, end, tmp)

# number of bases per bin, 3,031,053 bins
f2 <- fread(paste(wkdir, "bindata.1000.hg38.tsv.gz", sep="/")) %>% 
  filter(chr != "chrY") %>% 
  select(chr, start0, end, N, GC, k24gt90.sum, zooBasesFdr05, zooPrimate,
         cds, tssM50, els.sum:anchor.sum) 

# percentile rank of each bin
f3 <- f2 %>% 
  mutate_at(vars(N:anchor.sum), ~ percent_rank(.))

# interval join & drop outer edges
setDT(f1)
setkey(f2, chr, start0, end)
setkey(f3, chr, start0, end)

# output is fraction of bases per CNV
f4 <- foverlaps(f1, f2, type="any") %>% 
  select(chr, start0, end, i.start0, i.end, everything()) %>% 
  filter(between(start0, i.start0, i.end) & between(end, i.start0, i.end))
f4 <- f4 %>% 
  rename(GCtmp=GC, cdstmp=cds) %>% 
  group_by(tmp, i.start0, i.end) %>% 
  summarise(nBins = n(),
            Nbase = sum(N),
            GC = sum(GCtmp),
            k24 = sum(k24gt90.sum),
            zooConsAm = sum(zooBasesFdr05, na.rm = T),
            zooConsPr = sum(zooPrimate, na.rm = T),
            cds = sum(cdstmp),
            els = sum(els.sum),
            pls = sum(pls.sum),
            dhs = sum(dhs.sum),
            tfbs = sum(tfbs.sum),
            anchor = sum(anchor.sum)) %>% 
  ungroup %>% 
  mutate(across(Nbase:anchor, ~ .x/(i.end-i.start0))) %>% 
  select(-i.start0, -i.end)

write.table(f4, paste(wkdir, "UKBB_CNVs_for_AN_hg38.ALL.info", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")




#======== burden measured by length, CNV count, mammalian constraint, dosage sensitive genes ================================

library(dplyr)

wkdir="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden" ## path to conduct total CNV burden analyses in.

fam <- read.table(paste(wkdir, "UKBB_CNVs_for_AN_hg38.ALL.fam", sep="/"))
pheno <- read.table("/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat", header=TRUE)
fam$V5 <- factor(fam$V5)
fam$Array <- pheno$Array[match(fam$V1, pheno$cnvID)]
fam$Sex <- pheno$Sex[match(fam$V1, pheno$cnvID)]
fam$Age <- pheno$Age[match(fam$V1, pheno$cnvID)]
fam$V6 <- as.factor(fam$V6)
write.table(fam, paste(wkdir, "UKBB_CNVs_for_AN_hg38.ALL.pheno.fam", sep="/"), row.names=F, col.names=T, quote=F, sep="\t") ## save for future use.

# Function to calculate metrics
calculate_metrics <- function(cnv_data, fam, predictor_col, test_type) {
  
  print(predictor_col)
  m <- glm(as.formula(paste0("V6 ~ Sex + Age + Array + ", predictor_col)), data = fam, family = "binomial")
  or <- exp(cbind(OR = coef(m), confint(m))) %>% as.data.frame()
  
  avg <- fam %>% group_by(V6) %>% summarise_at(vars(count), list(mean = mean)) %>% as.data.frame()
  avg_case <- avg[2, 2]
  avg_cont <- avg[1, 2]
  
  vec <- c(test_type, predictor_col, 
           signif(or[5, 1], 3), 
           round(or[5, 2], 2), 
           round(or[5, 3], 2), 
           signif(summary(m)$coefficients[5, 4], 3), 
           signif(avg_case, 3), 
           signif(avg_cont, 3))
  
  return(vec)
}


# Function to process data for a specific test type
process_data <- function(wkdir, type, test_type, fam, predictor_col = "") {
  res <- data.frame()
  
  for (i in 1:length(type)) {
    cnv <- read.table(paste0(wkdir, "/UKBB_CNVs_for_AN_hg38.", type[i], ".cnv"), header = TRUE)
    
    if (test_type == "Length") {
      cnv$length <- as.numeric(cnv$BP2 - cnv$BP1)
      sum <- cnv %>% select(IID, length) %>% group_by(IID) %>% summarise(TotalLength=sum(length)) %>% as.data.frame()
      max <- signif(max(sum$TotalLength, na.rm=TRUE), 3)
      breaks <- seq(0, max, 100000)
      labels <- seq(1:(length(breaks)-1))
      sum$length_grouping <- cut(sum$TotalLength, seq(0, max, 100000), labels = labels)
      fam$length <- sum$length_grouping[match(fam$V1, sum$IID)]
      fam$length <- as.numeric(fam$length)
      fam$length <- ifelse(is.na(fam$length), 0, fam$length)
      predictor_col <- "length"
    } else if (test_type == "Count") {
      cnv_counts <- cnv %>% group_by(FID) %>% count() %>% as.data.frame()
      fam$count <- cnv_counts$n[match(fam$V1, cnv_counts$FID)]
      fam$count <- ifelse(is.na(fam$count), 0, fam$count)
      predictor_col <- "count"
    } else if (test_type == "Zoonomia") {
      cnv_info <- fread(paste(wkdir, "UKBB_CNVs_for_AN_hg38.", type[i], ".info", sep = "/"))
      cnv <- cnv %>% mutate(tmp = row_number())
      cnv_info <- merge(cnv, cnv_info, by = "tmp")
      predictor_col <- features[f]
    }
    
    vec <- calculate_metrics(cnv, fam, predictor_col, type[i])
    print(vec)
    res <- rbind(res, vec)
  }
  
  colnames(res) <- c("Type", "Test", "OR", "lower", "upper", "pval", "Avg_case", "Avg_cont")
  return(res)
}

type <- c("ALL", "DEL", "DUP")

# Process Length
res_length <- process_data(wkdir, type, "Length", fam)

# Process Count
res_count <- process_data(wkdir, type, "Count", fam)

# Process Zoonomia
features <- c("zooConsAm", "zooConsPr")
res_zoo <- process_data(wkdir, type, "Zoonomia", fam, features)

# Process dosage-sensitive protein-coding genes
res4 <- data.frame()
for (cnv_type in c("DEL", "DUP")) {
  cnv_sens <- read.table(paste(wkdir, paste("UKBB_CNVs_for_AN_hg38", cnv_type, "sensitive.cnv", sep = "."), sep = "/"), header = TRUE)
  fam$cnv <- as.factor(ifelse(fam$V1 %in% cnv_sens$FID, "yes", "no"))
  fam_cnvs <- fam[fam$cnv == "yes", ]
  
  m <- glm(V6 ~ Sex + Age + Array + cnv, data = fam, family = "binomial")
  or <- exp(cbind(OR = coef(m), confint(m))) %>% as.data.frame()
  OR <- signif(or[5, 1], 3)
  dp <- nchar(OR) - 2
  lower <- round(or[5, 2], 2)
  upper <- round(or[5, 3], 2)
  
  cnv_counts <- cnv_sens %>% group_by(FID) %>% count() %>% as.data.frame()
  fam$count <- cnv_counts$n[match(fam$V1, cnv_counts$FID)]
  fam$count <- ifelse(is.na(fam$count), 0, fam$count)
  fam$count <- as.numeric(fam$count)
  
  avg <- fam %>% group_by(V6) %>% summarise_at(vars(count), list(mean = mean)) %>% as.data.frame()
  avg_case <- avg[2, 2]
  avg_cont <- avg[1, 2]
  
  vec <- c(cnv_type, ifelse(cnv_type == "DEL", "Haploinsufficient", "Triplosensitive"), OR, lower, upper, 
           signif(summary(m)$coefficients[5, 4], 3), signif(avg_case, 3), signif(avg_cont, 3))
  res4 <- rbind(res4, vec)
}

colnames(res4) <- c("Type", "Test", "OR", "lower", "upper", "pval", "Avg_case", "Avg_cont")

# Combine all results
all_results <- rbind(res_length, res_count, res_zoo, res4)
all_results <- all_results[order(all_results$Type),]

# Write to CSV
write.table(all_results, paste(wkdir, "burden_results.csv", sep = "/"), row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)



#===============rare CNV burden across length and CNV type (all, del, dup)=====================================

fam <- read.table(paste(wkdir, "UKBB_CNVs_for_AN_hg38.ALL.pheno.fam", sep="/"))

lengths <- c("all", "first", "second", "third", "fourth")
svtype <- c("DEL", "DUP", "ALL")


## burden associations
res <- data.frame()
for (i in 1:length(lengths)){
  print(lengths[i])
  for (k in 1:length(svtype)){
    cnv <- read.table(paste0(wkdir, "/UKBB_CNVs_for_AN_hg38.", svtype[k],".", lengths[i], ".cnv"), header=T)
    fam$cnv <- as.factor(as.vector(ifelse(fam$V1%in%cnv$FID, "yes", "no")))
    m <- glm(V6 ~ Sex + Age + Array + cnv, data=fam, family="binomial")
    or <- exp(cbind(OR = coef(m), confint(m))) %>% as.data.frame()
    OR <- signif(or[5,1], 3)
    dp <- nchar(OR)-2
    lower <- round(or[5,2], 2)
    upper <- round(or[5,3], 2)
    pval <- signif(summary(m)$coefficients[5,4],3)
    cnv_cnt <- cnv %>% group_by(FID) %>% count() %>% as.data.frame()
    fam$count <- cnv_cnt$n[match(fam$V1, cnv_cnt$FID)]
    fam$count <- ifelse(is.na(fam$count), 0, fam$count)
    fam$count <- as.numeric(fam$count)
    avg <- fam %>% group_by(V6) %>% summarise_at(vars(count), list(mean=mean)) %>% as.data.frame()
    avg_case <- avg[2,2]
    avg_cont <- avg[1,2]
    vec <- c(svtype[k], lengths[i], OR, lower, upper, pval, signif(avg_case, 3), signif(avg_cont,3))
    print(vec)
    res <- rbind(res, vec)
  }
}   

    
#====== combine all burden results into one table 
colnames(res) <- c("Type", "Length","OR", "lower", "upper", "Pvalue", "Avg_Case", "Avg_Cont")

res$Length <- as.vector(ifelse(res$Length == "first", "20-100kB", ifelse(
  res$Length == "second", "100-200kB", ifelse(
    res$Length == "third", "200-500kB", ifelse(
      res$Length == "fourth", ">500kB", "All"
    )
  ))))
res$Type <- ifelse(res$Type=="ALL", "Deletions and Duplications", ifelse(res$Type=="DEL", "Deletions", "Duplications"))
res$Type <- factor(res$Type, levels=c("Deletions", "Duplications", "Deletions and Duplications"))
res2 <- res[with(res, order(Type)),]
res2$CI <- paste0("(", res2$lower, "-", res2$upper, ")")
write.table(res2, paste(wkdir, "UKBB_length_type.csv", sep="/"), row.names=FALSE, col.names=TRUE, sep =",", quote=TRUE)



#==== burden of rare CNVs split by CNV frequency (singleton up to 1%)==============================

freq <- fread(paste(wkdir, "/rare/UKBB_CNVs_for_AN_hg38.rare.freq.cnv", sep="/"))
freq <- as.data.frame(freq)
freq$V9 <- as.numeric(freq$V9)
freq$V3 <- as.factor(freq$V3)
summary(freq$V9)

freq$CNV_Count <- as.factor(as.vector(ifelse(freq$V9==1, "Singleton",
                                             ifelse(freq$V9>=2 &  freq$V9 <= 5, "2-5",
                                                    ifelse(freq$V9>=6 &  freq$V9 <= 10, "6-10",
                                                           ifelse(freq$V9>=11 &  freq$V9 <= 100, "11-100",
                                                                  ifelse(freq$V9>=101 &  freq$V9 <= 1000, "101-1000",
                                                                         ifelse(freq$V9>=1001 &  freq$V9 <= 1500, "1001-1500","1501-3872"))))))))


fam <- read.table(paste(wkdir, "UKBB_CNVs_for_AN_hg38.ALL.fam", sep="/"))
pheno <- read.table("/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat", header=TRUE)
fam$Array <- pheno$Array[match(fam$V1, pheno$V1)]
fam$Sex <- pheno$Sex[match(fam$V1, pheno$V1)]
fam$Age <- pheno$Age[match(fam$V1, pheno$V1)]

res <- data.frame()
unique <- levels(freq$CNV_Count)
for (u in 1:length(unique)){
  print(unique[u])
  cnvs <- freq %>% filter(CNV_Count==unique[u])
  fam$V5 <- factor(fam$V5)
  fam$V6 <- factor(fam$V6)
  fam$cnv <- as.factor(as.vector(ifelse(fam$V1%in%cnvs$V1, "yes", "no")))
  m <- glm(V6 ~ Sex + Age + Array + cnv, data = fam, family = "binomial")
  or <- exp(cbind(OR = coef(m), confint(m)))
  OR <- signif(or[5,1], 3)
  dp <- nchar(OR)-2
  lower <- round(or[5,2], 2)
  upper <- round(or[5,3], 2)
  pval <- signif(summary(m)$coefficients[5,4],3)
  cnv_cnt <- cnvs %>% group_by(V1) %>% count() %>% as.data.frame()
  fam$count <- cnv_cnt$n[match(fam$V1, cnv_cnt$V1)]
  fam$count <- ifelse(is.na(fam$count), 0, fam$count)
  fam$count <- as.numeric(fam$count)
  avg <- fam %>% group_by(V6) %>% summarise_at(vars(count), list(mean=mean)) %>% as.data.frame()
  avg_case <- avg[2,2]
  avg_cont <- avg[1,2]
  vec <- c(unique[u], OR, lower, upper, pval, avg_case, avg_cont)
  print(vec)
  res <- rbind(res, vec)
}

colnames(res) <- c("CNV_Count", "OR", "lower", "upper", "Pvalue", "Avg_Case", "Avg_Cont")
write.table(res, paste(wkdir, "rare/UKB_rare_counts.csv", sep="/"), col.names = T, row.names = F, sep=",", quote=F)


