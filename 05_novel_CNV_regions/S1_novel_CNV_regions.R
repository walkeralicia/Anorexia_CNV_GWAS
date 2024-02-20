
#========== This script identifies novel AN-associated disease-risk CNV regions and intersecting genes =========================================================

#== Load required R libraries
library(dplyr)
library(data.table)

#== Set up file paths ===
wkdir="/Anorexia/UKB/burden_analysis/SAIGE"  ## path to SAIGE result directory
data_path="data" ## data path to cytoBand.txt, and geneMatrix.tsv.

#======= Read in Genomic cytoband information===============================================
cytoband <- fread(paste(data_path, "cytoBand.txt", sep="/"), header=FALSE) %>% as.data.frame()
cytoband$CHR <- as.numeric(gsub("chr", "", cytoband$V1))
cytoband$CHR <- gsub("X", "23", cytoband$CHR)
cytoband$CHR <- as.numeric(cytoband$CHR)

#======= Read in information for all SNPs genotyped within UKB individuals used for CNV calling ===========================
bim <- fread("UKB.bim", header = F)
bim <- as.data.frame(bim)

#=== DELETIONS ==========================================================================================================================

#=== Extract all UKB CNV deletions (both rare and common CNVs) ===========
system("cfile=/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38")
system(paste("plink --noweb --cfile ${cfile} --cnv-del --cnv-write --out ${cfile}.DEL"))
system(paste("plink --noweb --cfile ${cfile}.DEL --cnv-make-map --out ${cfile}.DEL"))

#=== Identify all unique CNV deletion breakpoints (rare and common) ====================
del_cnvs <- fread("/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.DEL.cnv")
breakpoints <- unique(c(del_cnvs$BP1, del_cnvs$BP2))

#=== Load SAIGE rare CNV breakpoint GWAS results for deletions ====================
dels <- read.table(paste(wkdir, "step2/DEL_SAIGE_step2_AN_rare.txt", sep="/"), header=TRUE)
colnames(dels)[3] <- "ID"

#=== Identify novel AN-association CNV regions looping by chromosome================
cnvr_dels_3000kb <- data.frame()

for (c in 1:23){
  print(c)
  dels_chr <- dels %>% dplyr::filter(CHR==c)
  sig_c <- dels_chr %>% filter(p.value<0.05)
  sig_c_ord <- sig_c[order(sig_c$p.value, decreasing=F),]
  sig_c_ord <- sig_c_ord %>% filter(ID%in%breakpoints)
  sic_c_ord <- sig_c_ord %>% filter(ID%in%dels$ID)
  sig_rs <- sig_c_ord$ID
  bim_c <- bim[bim$V1==c,]
  if (length(sig_rs)==0) next
  
  ## loop over sig SNPs that are not already in a CNVR
  for (i in 1:nrow(sig_c_ord)){
    rs <- sig_c_ord[i,3]
    chr <- sig_c_ord[i,1]
    print(paste0("Starting rs ",rs," chr", chr, " number", i))
    
    if (!sig_c_ord[i,3]%in%sig_rs | length(sig_rs)==0) next
    print("rs not in CNVR yet")
    
    ## save temp output file
    rs_file <- paste0(wkdir, "/CNVR/DEL/Chr", chr, "_rs", rs)
    fwrite(data.frame(ID = rs), rs_file, col.names = F, row.names = F, quote = F, sep = "\t")
    
    # Run PLINK: detect probes in LD
    input <- paste0("/Anorexia/UKB/plink_files/bfiles/cnv_to_PLINK_DEL_chr", chr)
    system(paste("plink",
                 "--bfile", input,
                 "--show-tags", rs_file,
                 "--tag-kb 3000 --tag-r2 0.5",
                 "--out", rs_file))
    unlink(rs_file); unlink(paste0(rs_file, ".log"))
    
    # Detect the boundaries of the new CNVR
    cnvr <- as.data.frame(fread(paste0(rs_file, ".tags"), header = F, col.names = "ID"))
    cnvr <- left_join(cnvr, dels_chr, by = "ID")
    cnvr <- cnvr %>% filter(ID%in%breakpoints)
    cnvr <- cnvr %>% filter(ID%in%dels$ID)
    rs_pval <- signif(cnvr[which(cnvr$ID==rs), 13],2)
    start <- min(cnvr$ID)
    end <- max(cnvr$ID)
    length <- round((end-start)/1000, 1)
    region <- paste0("chr", chr, ":", start, "-",end)
    
    # Extract rare CNVb burden results
    num <- nrow(cnvr)
    bim_c_cnv <- bim_c %>% filter(V4 >= start & V4 <= end)
    num_probes <- nrow(bim_c_cnv)
    cnvr$freq <- (100*cnvr$AC_Allele2)/387190 ## change to the sample size
    freq <- signif(cnvr[which(cnvr$ID==rs), 24],2)
    beta <- cnvr[which(cnvr$ID==rs),9]
    OR <- signif(exp(beta),2)
    dp <- nchar(OR)-2
    se <- cnvr[which(cnvr$ID==rs),10]
    LCI <- round(exp(beta-1.96*se), dp)
    UCI <- round(exp(beta+1.96*se), dp)
    CI <- paste0("(", LCI, "-", UCI, ")")
    n_case <- cnvr[which(cnvr$ID==rs), 21]
    n_cont <- cnvr[which(cnvr$ID==rs), 23]
    cytoband_chr = cytoband[cytoband$CHR==chr,]
    cytoband_chr$pos = ifelse(cytoband_chr$V2 <= rs & cytoband_chr$V3 >= rs, "yes", "no")
    cyto <- cytoband_chr[which(cytoband_chr$pos=="yes"), "V4"]
    cyto <- paste0(chr, cyto)
    vec <- c(cyto, region, "Deletion",length, num_probes, num, rs, rs_pval, OR, CI, freq, n_case, n_cont)
    names(vec) <- c("Cytoband","CNVR","Genotype", "Length", "Num_Probes", "Num_Breakpoints", "Top_Breakpoint", "Pval", "OR", "CI", "CNV_Freq", "N_Case", "N_Cont")
    print(vec)
    cnvr_dels_3000kb <- rbind(cnvr_dels_3000kb, vec)
    sig_rs <- setdiff(sig_rs, cnvr$ID)
  }
}

colnames(cnvr_dels_3000kb) <- c("Cytoband", "CNVR", "Genotype", "Length", "Num_Probes","Num_Breakpoints", "Top_Breakpoint", "Pval", "OR", "CI", "CNV_Freq", "N_Case", "N_Cont")
dels_ord <- cnvr_dels_3000kb[order(cnvr_dels_3000kb$Pval),]
dels_ord$OR <- as.numeric(dels_ord$OR)
dels_ord$Num_Probes <- as.numeric(dels_ord$Num_Probes)
dels_ord$Length <- as.numeric(dels_ord$Length)
## Only keep CNVR longer than 20kb and that have at least 10 probes
dels_ord_filter <- dels_ord %>% filter(Num_Probes>=10 & Length >=20)
dels_final <- dels_ord_filter


#=== DUPLICATIONS ==========================================================================================================================

#=== Extract all UKB CNV duplications (rare and common CNVs) ===========
system("cfile=/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38")
system(paste("plink --noweb --cfile ${cfile} --cnv-dup --cnv-write --out ${cfile}.DUP"))
system(paste("plink --noweb --cfile ${cfile}.DUP --cnv-make-map --out ${cfile}.DUP"))

#=== Identify all unique CNV deletion breakpoints (rare and common) ====================
dup_cnvs <- fread("/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.DUP.cnv")
breakpoints <- unique(c(dup_cnvs$BP1, dup_cnvs$BP2))

#=== Load SAIGE rare CNV breakpoint GWAS results for deletions ====================
dups <- read.table(paste(wkdir, "step2/DUP_SAIGE_step2_AN_rare.txt", sep="/"), header=TRUE)
colnames(dups)[3] <- "ID"

#=== Identify novel AN-association CNV regions looping by chromosome================
cnvr_dups_3000kb <- data.frame()

## loop by chromosomes
for (c in 1:23){
  dels_chr <- dels %>% dplyr::filter(CHR==c)
  sig_c <- dels_chr %>% filter(p.value<0.05)
  sig_c_ord <- sig_c[order(sig_c$p.value, decreasing=F),]
  sig_c_ord <- sig_c_ord %>% filter(ID%in%breakpoints)
  sic_c_ord <- sig_c_ord %>% filter(ID%in%dups$ID)
  sig_rs <- sig_c_ord$ID
  if (length(sig_rs)==0) next
  
  ## loop over sig SNPs that are not already in a CNVR
  for (i in 1:nrow(sig_c_ord)){
    rs <- sig_c_ord[i,3]
    chr <- sig_c_ord[i,1]
    print(paste0("Starting rs ",rs," chr", chr, " number", i))
    
    if (!sig_c_ord[i,3]%in%sig_rs | length(sig_rs)==0) next
    print("rs not in CNVR yet")
    
    ## save temp output file
    rs_file <- paste0(wkdir, "/CNVR/DUP/Chr", chr, "_rs", rs)
    fwrite(data.frame(ID = rs), rs_file, col.names = F, row.names = F, quote = F, sep = "\t")
    
    # Run PLINK: detect probes in LD
    input <- paste0("/Anorexia/UKB/plink_files/bfiles/cnv_to_PLINK_DUP_chr", chr)
    system(paste("plink",
                 "--bfile", input,
                 "--show-tags", rs_file,
                 "--tag-kb 3000 --tag-r2 0.5",
                 "--out", rs_file))
    unlink(rs_file); unlink(paste0(rs_file, ".log"))
    
    # Detect the boundaries of the CNVR
    cnvr <- as.data.frame(fread(paste0(rs_file, ".tags"), header = F, col.names = "ID"))
    cnvr <- left_join(cnvr, dups_chr, by = "ID")
    cnvr <- cnvr %>% filter(ID%in%breakpoints)
    cnvr <- cnvr %>% filter(ID%in%dups$ID)
    rs_pval <- signif(cnvr[which(cnvr$ID==rs), 13],2)
    start <- min(cnvr$ID)
    end <- max(cnvr$ID)
    length <- round((end-start)/1000, 1)
    region <- paste0("chr", chr, ":", start, "-",end)
    num <- nrow(cnvr)
    bim_c_cnv <- bim_c %>% filter(V4 >= start & V4 <= end)
    num_probes <- nrow(bim_c_cnv)
    cnvr$freq <- (100*cnvr$AC_Allele2)/387190 ## change to the sample size
    freq <- signif(cnvr[which(cnvr$ID==rs), 24],2)
    beta <- cnvr[which(cnvr$ID==rs),9]
    OR <- signif(exp(beta),2)
    dp <- nchar(OR)-2
    se <- cnvr[which(cnvr$ID==rs),10]
    LCI <- round(exp(beta-1.96*se), dp)
    UCI <- round(exp(beta+1.96*se), dp)
    CI <- paste0("(", LCI, "-", UCI, ")")
    n_case <- cnvr[which(cnvr$ID==rs), 21]
    n_cont <- cnvr[which(cnvr$ID==rs), 23]
    cytoband_chr = cytoband[cytoband$CHR==chr,]
    cytoband_chr$pos = ifelse(cytoband_chr$V2 <= rs & cytoband_chr$V3 >= rs, "yes", "no")
    cyto <- cytoband_chr[which(cytoband_chr$pos=="yes"), "V4"]
    cyto <- paste0(chr, cyto)
    vec <- c(cyto, region, "Duplication", length, num_probes, num, rs, rs_pval, OR, CI, freq, n_case, n_cont)
    names(vec) <- c("Cytoband", "CNVR","Genotype", "Length", "Num_Probes", "Num_Breakpoints", "Top_Breakpoint", "Pval", "OR", "CI", "CNV_Freq", "N_Case", "N_Cont")
    print(vec)
    cnvr_dups_3000kb <- rbind(cnvr_dups_3000kb, vec)
    sig_rs <- setdiff(sig_rs, cnvr$ID)
  }
}

colnames(cnvr_dups_3000kb) <- c("Cytoband", "CNVR", "Genotype", "Length", "Num_Probes","Num_Breakpoints", "Top_Breakpoint", "Pval", "OR", "CI", "CNV_Freq", "N_Case", "N_Cont")
dups_ord <- cnvr_dups_3000kb[order(cnvr_dups_3000kb$Pval),]
dups_ord$OR <- as.numeric(dups_ord$OR)
dups_ord$Num_Probes <- as.numeric(dups_ord$Num_Probes)
dups_ord$Length <- as.numeric(dups_ord$Length)
## Only keep CNVR longer than 20kb and that have at least 10 probes
dups_ord_filter <- dups_ord %>% filter(Num_Probes>=10 & Length >=20)
final_dup <- dups_ord_filter

#====== Combine CNVR deletions and duplications into one dataset=============================
both <- rbind(dels_final, final_dup)
both <- unique(both)
both_risk <- both %>% filter(OR>1 & N_Case > 0) 
both_risk_ord <- both_risk[order(both_risk$Pval),]

#======= Load gene annotation file ============================================================
genes <- read.table(paste0(data_path, "/geneMatrix.tsv"), header=TRUE, fill=TRUE)
genes <- genes %>% filter(gene_type=="protein_coding")
genes3 <- genes[, c(8,9,10,1)]
fwrite(genes3, paste("CNVR/Genes",sep="/"), col.names=F, sep="\t")
genes2 <- genes[, c(8,9,10,2)]
fwrite(genes2, paste(wkdir, "CNVR/Genes2", sep="/"), col.names=F, sep="\t")

#===== Intersect novel CNVRs with genes ==========================================================
cnvs <- both_risk_ord
cnvs$chrom <- sapply(strsplit(cnvs$CNVR, ":"), "[[", 1)
cnvs$Start <- as.numeric(sapply(strsplit(sapply(strsplit(cnvs$CNVR, ":"), "[[", 2), "-"), "[[", 1))
cnvs$End <- as.numeric(sapply(strsplit(sapply(strsplit(cnvs$CNVR, ":"), "[[", 2), "-"), "[[", 2))
cnvs2 <- cnvs %>% select(chrom, Start, End, CNVR)

both_risk_ord$Genes <- NA
both_risk_ord$Gene_Names <- NA
for (i in 1:nrow(both_risk_ord)){
  cnvr <- cnvs2[i,]
  fwrite(cnvr, "TMP1", col.names = F, sep = "\t")
  system(paste("bedtools intersect -wao -a TMP1",
               "-b /Anorexia/UKB/burden_analysis/SAIGE/CNVR/Genes",
               "> TMP2"))
  g <- fread("TMP2")
  all <- paste(g$V8,collapse=" ")
  both_risk_ord$Genes[i] <- all
  
  system(paste("bedtools intersect -wao -a TMP1",
               "-b /Anorexia/UKB/burden_analysis/SAIGE/CNVR/Genes2",
               "> TMP2"))
  g <- fread("TMP2")
  all <- paste(g$V8,collapse=" ")
  both_risk_ord$Gene_Names[i] <- all
}

both_risk_ord$Genes <- ifelse(both_risk_ord$Genes == ".", NA, both_risk_ord$Genes)
both_risk_ord$Gene_Names <- ifelse(both_risk_ord$Gene_Names == ".", NA, both_risk_ord$Gene_Names)

#==== Save all results==========================
write.table(both_risk_ord, paste(wkdir, "CNVR/all_CNVR_risk_genes.csv", sep="/"), col.names=T, row.names=F, quote=F, sep="," )





