

#================= This script evaluates total rare CNV burden ====================================================

#============= Load CNV burden function, file paths, and read in phenotypic data for covariates ========================

# Load required R libraries
library(dplyr)
library(data.table)

# File paths and burden function
wkdir<-"/Anorexia/UKB/burden_analysis/overall_burden" ## path to conduct total CNV burden analyses in.
cnv_path<-"/Anorexia/UKB/rare_cnvs" ## path to rare CNV files
pheno_path<-"/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat" ## path to phenotype file for covariates
cfile<-"UKBB_CNVs_for_AN_hg38"
source("calculate_burden_function.R") ## Burden function


#==========CNV burden measured by CNV length, count, prop. of mammalian constraint bases, and dosage sensitive genes =================

pheno <- read.table(pheno_path, header=T)

# Function to process data for a specific test type
process_data <- function(wkdir, cnv_type, test_type, pheno) {
  res <- list()
  
  for (i in 1:length(cnv_type)) {
    cnv <- read.table(paste0(cnv_path, "/", cfile,".", cnv_type[i], ".cnv"), header = TRUE)
    
    if (test_type == "Length") {
      # CNV length burden
      cnv$length <- as.numeric(cnv$BP2 - cnv$BP1)
      sum <- cnv %>% select(IID, length) %>% group_by(IID) %>% summarise(TotalLength=sum(length)) %>% as.data.frame()
      max <- signif(max(sum$TotalLength, na.rm=TRUE), 3)
      breaks <- seq(0, max, 100000)
      labels <- seq(1:(length(breaks)-1))
      sum$length_grouping <- cut(sum$TotalLength, seq(0, max, 100000), labels = labels)
      pheno$Length <- sum$length_grouping[match(pheno$V1, sum$IID)]
      pheno$Length <- as.numeric(pheno$Length)
      pheno$Length <- ifelse(is.na(pheno$Length), 0, pheno$Length)
      
    } else if (test_type == "Count") {
      # CNV count burden 
      cnv_counts <- cnv %>% group_by(FID) %>% count() %>% as.data.frame()
      pheno$Count <- cnv_counts$n[match(pheno$V1, cnv_counts$FID)]
      pheno$Count <- ifelse(is.na(pheno$Count), 0, pheno$Count)
      
    } else if (test_type == "zooConsAm") {
      # All mammalian constraint
      cnv_info <- fread(paste0(cnv_path, "/", cfile, ".ALL.info"))
      cnv <- cnv %>% mutate(tmp = row_number())
      cnv_info <- merge(cnv, cnv_info, by = "tmp")
      avg <- cnv_info %>% select(IID, zooConsAm) %>% group_by(IID) %>% summarise(avg=mean(zooConsAm)) %>% as.data.frame()
      pheno$zooConsAm <- avg$avg[match(pheno$V1, avg$IID)]
      pheno$zooConsAm <- ifelse(is.na(pheno$zooConsAm), 0, pheno$zooConsAm)
      
    } else if (test_type == "zooConsPr") {
      # Primate constraint
      cnv_info <- fread(paste0(cnv_path, "/", cfile, ".ALL.info"))
      cnv <- cnv %>% mutate(tmp = row_number())
      cnv_info <- merge(cnv, cnv_info, by = "tmp")
      avg <- cnv_info %>% select(IID, zooConsPr) %>% group_by(IID) %>% summarise(avg=mean(zooConsPr)) %>% as.data.frame()
      pheno$zooConsPr <- avg$avg[match(pheno$V1, avg$IID)]
      pheno$zooConsPr <- ifelse(is.na(pheno$zooConsPr), 0, pheno$zooConsPr)
      
    } else if (test_type == "dosage_sensitive") {
      # Dosage-sensitive constraint
      cnv <- read.table(paste0(cnv_path, "/", cfile, ".", cnv_type[i], ".sensitive.cnv"), header = TRUE)
      pheno$dosage_sensitive <- ifelse(pheno$V1 %in% cnv$FID, 1, 0)
      
    }
    vec <- calculate_burden(cnv, pheno, test_type, cnv_type[i])
    res[[i]] <- vec
  }
  res_all <- do.call(rbind, res)
  return(res_all)
}

cnv_type <- c("ALL", "DEL", "DUP")

# Process tests
res_length <- process_data(wkdir, cnv_type, "Length", pheno)
res_count <- process_data(wkdir, cnv_type, "Count", pheno)
res_zooConsAm <- process_data(wkdir, cnv_type, "zooConsAm", pheno)
res_zooConsPr <- process_data(wkdir, cnv_type, "zooConsPr", pheno)
res_ds <- process_data(wkdir, cnv_type[2:3], "dosage_sensitive", pheno)

# Combine all results
all_results <- rbind(res_length, res_count, res_zooConsAm, res_zooConsPr, res_ds) %>% as.data.frame()
all_results <- all_results[order(all_results$Type),]

# Write to CSV
write.table(all_results, paste(wkdir, "UKBB_overall_burden_results.csv", sep = "/"), row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)



#===============rare CNV burden across length and CNV type (all, del, dup)=====================================

lengths <- c("all", "first", "second", "third", "fourth")
svtype <- c("DEL", "DUP", "ALL")

## burden associations
length_list <- list()
for (i in 1:length(lengths)){
  print(lengths[i])
  svtype_list <- list()
  for (k in 1:length(svtype)){
    print(svtype[k])
    cnv <- read.table(paste0(cnv_path, "/", cfile, ".", svtype[k],".", lengths[i], ".cnv"), header=T)
    pheno <- read.table(pheno_path, header=T)
    pheno$cnv <- as.vector(ifelse(pheno$V1%in%cnv$FID, 1, 0))
    colnames(pheno)[7] <- lengths[i]
    predictor_col =  lengths[i]
    vec <- calculate_burden(cnv, pheno, lengths[i], svtype[k])
    svtype_list[[k]] <- vec
  }
  svtype_list_all <- do.call(rbind, svtype_list)
  length_list[[i]] <- svtype_list_all
}

length_list_all <- do.call(rbind, length_list) %>% as.data.frame()


# combine all burden results into one table 

res <- length_list_all
res$Test <- as.vector(ifelse(res$Test == "first", "20-100kB", ifelse(
  res$Test == "second", "100-200kB", ifelse(
    res$Test == "third", "200-500kB", ifelse(
      res$Test == "fourth", ">500kB", "All"
    )
  ))))
res$Type <- ifelse(res$Type=="ALL", "Deletions and Duplications", ifelse(res$Type=="DEL", "Deletions", "Duplications"))
res$Type <- factor(res$Type, levels=c("Deletions", "Duplications", "Deletions and Duplications"))
res2 <- res[with(res, order(Type)),]
colnames(res2)[2] <- "Length"
write.table(res2, paste(wkdir, "UKBB_length_type_results.csv", sep="/"), row.names=FALSE, col.names=TRUE, sep =",", quote=TRUE)



#========== rare CNV burden split by CNV frequency (singleton up to 1%)==============================

freq <- fread(paste(cnv_path, "/", cfile, ".rare.freq.cnv", sep="/"))
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





unique <- c("Singleton", "2-5", "6-10", "11-100", "101-1000", "1001-1500", "1501-3872")
count_list <- list()

for (u in 1:length(unique)){
  print(unique[u])
  pheno <- read.table(pheno_path, header=T)
  cnv <- freq %>% filter(CNV_Count==unique[u])
  pheno$freq <- as.vector(ifelse(pheno$V1%in%cnv$V1, 1, 0))
  vec <- calculate_burden(cnv, pheno, "freq", "ALL")
  count_list[[u]] <- vec
}
count_list_all <- do.call(rbind, count_list) %>% as.data.frame()
count_list_all$Test <- unique

write.table(count_list_all, paste(wkdir, "UKBB_rare_counts.csv", sep="/"), col.names = T, row.names = F, sep=",", quote=F)


