

#================ This script extracts UKB CNVs intersecting novel ANGI CNV regions associated with AN status =====================================

#== Load R libraries ==
library(dplyr)
library(data.table)

#== Set up paths ==
wkdir <- "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE/step2"
ANGI_data<-"ANGI_data" ## path to folder with ANGI meta-analyses results

#== Load data ==
cnvrs <- read.csv(paste(ANGI_data, "ANGI_CNVR_zoonomia.csv", sep ="/"), header=T)
saige_dels <- read.table(paste(wkdir, "DEL_SAIGE_step2_AN_rare.txt", sep="/"), header=T)
saige_dups <- read.table(paste(wkdir, "DUP_SAIGE_step2_AN_rare.txt", sep="/"), header=T)

# begin UKB CNV overlap with novel ANGI CNV regions ================================================================

cnvs <- data.frame()
for (i in 1:nrow(cnvrs)){
  
  ## extract ANGI CNVR boundaries ======================================================
  geno <- cnvrs[i, "Genotype"]
  CNVR <- cnvrs[i, "CNVR"]
  c <- gsub("chr", "", strsplit(CNVR, ":")[[1]][1])
  start <- as.numeric(strsplit(strsplit(CNVR, ":")[[1]][2], "-")[[1]][1])
  end <- as.numeric(strsplit(strsplit(CNVR, ":")[[1]][2], "-")[[1]][2])
  
  ## extract UKB SAIGE CNVb GWAS results within CNVR with the corresponding genotype ====
  if (geno=="Deletion"){
    cnvr_saige <- saige_dels
  } else {
    cnvr_saige <- saige_dups
  }
  cnvr_saige <- cnvr_saige %>% filter(CHR==c) %>% filter(MarkerID >= start & MarkerID <= end)
  num_breakpoints <- nrow(cnvr_saige)
  
  ## extract UKB CNVs that intersect the CNVR ===================================================
  cnvr_tmp <- c(c, start, end, "tmp")
  cnvr_tmp <- data.frame(t(unlist(cnvr_tmp)))
  fwrite(cnvr_tmp, "TMP1", col.names = F, sep = "\t")
  if (geno=="Deletion"){
    system(paste("plink --noweb --cfile UKBB_CNVs_for_AN_hg38.DEL --cnv-intersect TMP1 --cnv-write --out TMP2"))
  } else {
    system(paste("plink --noweb --cfile UKBB_CNVs_for_AN_hg38.DUP --cnv-intersect TMP1 --cnv-write --out TMP2"))
  }
  tracks <- fread("TMP2.cnv")
  fam <- fread("TMP2.fam")
  tracks$status <- fam$V6[match(tracks$FID, fam$V1)]
  tracks_controls <- tracks %>% filter(status==1)
  num_cont <- nrow(tracks_controls)
  num_case <- nrow(tracks)-nrow(tracks_controls)
  
  ## extract UKB SAIGE result for the top ANGI breakpoint within corresponding CNVR ============
  cnvr_top <- cnvr_saige %>% filter(p.value ==min(cnvr_saige$p.value))
  cnvr_top <- cnvr_top[1,]
  top_bp <- cnvr_top[, "MarkerID"]
  top_pval <- signif(cnvr_top[, "p.value"],2)
  beta <- cnvr_top[, "BETA"]
  OR <- signif(exp(beta), 2)
  dp <- nchar(OR)-2
  se <- cnvr_top[, "SE"]
  LCI <- round(exp(beta-1.96*se), dp)
  UCI <- round(exp(beta+1.96*se), dp)
  CI <- paste0("(", LCI, "-", UCI, ")")
  
  ## write table of results ==============================================
  cnvrs_angi <- cnvrs %>% select(Cytoband, CNVR, Genotype, Length, Num_Probes, Num_Breakpoints,
                                 N_Case_CNV, N_Cont_CNV, 
                                 Top_Breakpoint, Pval, OR, CI)
  cnvrs_constraint <- cnvrs %>% select(Genes, Gene_Names, phyloPConst_prop, Mean_phyloP,Max_phyloP,
                                       phastConst_prop, Mean_phastCons, Max_phastCons)
  cnvrs_ukb <- data.frame(t(c(num_breakpoints, num_case, num_cont, top_bp, top_pval, OR, CI)))
  names(cnvrs_ukb) <- c("UKB_Num_Breakpoints", "UKB_Num_Case", "UKB_Num_Cont", "UKB_Top_BP",
                        "UKB_Pval", "UKB_OR", "UKB_CI")
  res <- cbind(cnvrs_angi[i,],cnvrs_ukb,cnvrs_constraint[i,])
  cnvs <- rbind(cnvs,res)
  
}

write.csv(cnvs, paste(wkdir, "angi_ukb_cnvrs.csv", sep="/"))
