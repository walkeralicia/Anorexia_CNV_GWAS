
#============== This script prepares the pleiotropic CNV list ============================

#== Load required R libraries =======
library(dplyr)
library(readxl)
library(stringr)

#== Set up paths ====

data_path <- "data" ## path to 1-s2.0-S0092867422007887-mmc4.xlsx
wkdir <- "/Anorexia/UKB/burden_analysis/drCNVs" ## path to conduct disease-risk CNV burden analyses in

#== Format CNV list into .bed format ===========

collins <- read_excel(paste(data_path, "1-s2.0-S0092867422007887-mmc4.xlsx", sep="/"), sheet=1)
collins <- as.data.frame(collins)
dels <- collins %>% filter(`CNV Type`=="DEL")
dels <- dels %>% select(Chrom, Start, End, `Segment ID`)
write.table(dels, paste(wkdir, "collins_dels.txt", sep="/"), sep="\t", col.names=F, row.names=F, quote=F)
dups <- collins %>% filter(`CNV Type`=="DUP")
dups <- dups %>% select(Chrom, Start, End, `Segment ID`)
write.table(dups, paste(wkdir, "collins_dups.txt", sep="/"), sep="\t", col.names=F, row.names=F, quote=F)