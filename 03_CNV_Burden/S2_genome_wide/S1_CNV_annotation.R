

#======== This script annotates rare CNVs to dosage-sensitive genes and to mammalian constraint measures ===============================

#== load R packages ===
library(dplyr)
library(tidyverse)
library(data.table)
library(skimr)
library(rstatix)
library(stringr)
library(ggplot2)

#== Set up paths =====

data_path<-"data" ## path to data folder with bindata.1000.hg38.tsv.gz, geneMatrix.tsv, and collins2022-TableS7-pHaplo-pTriplo.tsv files.
cfile_path<-"/Anorexia/UKB/rare_cnvs" ## path to hg38 UKB rare CNVs
file_name<-"UKBB_CNVs_for_AN_hg38" ## cfile name
p<-"plink" ## plink software
drCNVS<-"/Anorexia/UKB/burden_analysis/drCNVs" ## path to conduct disease-risk CNV burden analyses in


#======== Extract collins et al. dosage sensitive genes ===========================


# Reformat list of genes, split by CNV type (i.e., triplosensitive (dups) and haploinsufficient (dels))
collins <- read.table(paste(data_path,"collins2022-TableS7-pHaplo-pTriplo.tsv", sep="/"), header=TRUE)
genes <- read.table(paste(data_path, "geneMatrix.tsv", sep="/"), header=TRUE, fill=TRUE)
genes2 <- genes[, c("ensgid", "hg38h0", "h1","h2")]
both <- merge(genes2, collins, by="ensgid", all.y=TRUE)
haplo <- both %>% filter(pHaplo >= 0.86)
triplo <- both %>% filter(pTriplo >= 0.96)
haplo2 <- haplo %>% select(hg38h0, h1, h2, ensgid)
haplo2$hg38h0 <- gsub("chr", "", haplo2$hg38h0)
triplo2 <- triplo %>% select(hg38h0, h1, h2, ensgid)
triplo2$hg38h0 <- gsub("chr", "", triplo2$hg38h0)
write.table(triplo2, paste(drCNVs,"triplosensitive_list.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
write.table(haplo2, paste(drCNVs,"haploinsufficient_list.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")


#============ Annotate rare UKB CNVs with disease-risk Collins dosage-sensitive genes ==================================

# Intersections
system(paste(p, "--noweb --cfile", paste(cfile_path, paste0(file_name, ".DUP"), sep = "/"), "--cnv-intersect", paste(drCNVS, "triplosensitive_list.txt", sep = "/"), "--cnv-write", "--out", paste(cfile_path, paste0(file_name, ".DUP.sensitive"), sep = "/")))

system(paste(p, "--noweb --cfile", paste(cfile_path, paste0(file_name, ".DUP.sensitive"), sep = "/"), "--cnv-make-map", "--out", paste(cfile_path, paste0(file_name, ".DUP.sensitive"), sep = "/")))

system(paste(p, "--noweb --cfile", paste(cfile_path, paste0(file_name, ".DEL"), sep = "/"), "--cnv-intersect", paste(drCNVS, "haploinsufficient_list.txt", sep = "/"), "--cnv-write", "--out", paste(cfile_path, paste0(file_name, ".DEL.sensitive"), sep = "/")))

system(paste(p, "--noweb --cfile", paste(cfile_path, paste0(file_name, ".DEL.sensitive"), sep = "/"), "--cnv-make-map", "--out", paste(cfile_path, paste0(file_name, ".DEL.sensitive"), sep = "/")))


#================ Annotate rare CNVs with mammalian constraint (1K bins from Pat Sullivan)======================

a <- fread(paste0(cfile_path, "/", file_name, ".ALL.cnv")) %>% 
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
f2 <- fread(paste(data_path, "bindata.1000.hg38.tsv.gz", sep="/")) %>% 
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

write.table(f4, paste0(cfile_path, "/", file_name, ".ALL.info"), row.names=F, col.names=T, quote=F, sep="\t")

