


#======================= This script prepares developmental CNV list ============================================================================

# Load required R libraries
library(dplyr)
library(readxl)
library(stringr)
library(data.table)
options(scipen=999)

# Format CNV list into .bed format and separate by CNV type (duplications and deletions)
cnv_path <- "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs"
cnvs <- read_excel(paste(cnv_path, "CNV_List.xlsx", sep="/"), sheet=1)
cnvs <- as.data.frame(cnvs)
cnvs$ID <- 1:nrow(cnvs)
cnvs$chrom <- sapply(strsplit(cnvs$Location, ":"), "[[", 1)
cnvs$chrom <- ifelse(cnvs$chrom == "X", 23, cnvs$chrom)
cnvs$Start <- sapply(strsplit(sapply(strsplit(cnvs$Location, ":"), "[[", 2), "-"), "[[", 1)
cnvs$End <- sapply(strsplit(sapply(strsplit(cnvs$Location, ":"), "[[", 2), "-"), "[[", 2)
dels <- cnvs %>% filter(Genotype=="Deletion")
dels <- dels %>% select(chrom, Start, End, ID)
write.table(dels, paste(cnv_path, "devCNVs_dels.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
dups <- cnvs %>% filter(Genotype=="Duplication")
dups <- dups %>% select(chrom, Start, End, ID)
write.table(dups, paste(cnv_path, "devCNVs_dups.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


# Format NRXN1 exon boundaries into .bed format 
df <- fread(paste(cnv_path, "NRXN1_exons", sep="/"),fill=TRUE, sep=",") %>% as.data.frame()
exons <- data.frame()
for (i in 1:nrow(df)){
  chr=2
  start=strsplit(df[1,"exonStarts"], split=",")[[1]][i]
  end=strsplit(df[1,"exonEnds"], split=",")[[1]][i]
  index=paste0(61, "_", i)
  vec <- c(chr, start, end, index)
  exons <- rbind(exons, vec)
}
write.table(exons, paste(cnv_path, "NRXN1_exon_boundaries.bed", sep="/"), col.names=F, row.names=F, sep="\t", quote=F)
