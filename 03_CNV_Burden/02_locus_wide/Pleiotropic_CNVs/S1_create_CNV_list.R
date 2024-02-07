
#============== this script prepares the pleiotropic CNV list ============================

# load required R libraries =======
library(dplyr)
library(readxl)
library(stringr)
options(scipen=999)

# format CNV list into .bed format ===========
wkdir <- "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs"
collins <- read_excel(paste(wkdir, "1-s2.0-S0092867422007887-mmc4.xlsx", sep="/"), sheet=1)
collins <- as.data.frame(collins)
dels <- collins %>% filter(`CNV Type`=="DEL")
dels <- dels %>% select(Chrom, Start, End, `Segment ID`)
write.table(dels, paste(wkidr, "collins_dels.txt", sep="/"), sep="\t", col.names=F, row.names=F, quote=F)
dups <- collins %>% filter(`CNV Type`=="DUP")
dups <- dups %>% select(Chrom, Start, End, `Segment ID`)
write.table(dups, paste(wkdir, "collins_dups.txt", sep="/"), sep="\t", col.names=F, row.names=F, quote=F)