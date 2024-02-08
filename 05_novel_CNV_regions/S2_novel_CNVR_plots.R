


#==== This script creates plots of the novel AN-risk CNVRs with intersecting individual-level CNV tracks and gene tracks underneath ================

# Load R libraries ==============
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
wkdir="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE"

# Read in AN-risk CNVRs =====================
cnvr <- read.csv(paste(wkdir, "CNVR/all_CNVR_risk_genes.csv", sep="/"), header=T)
cnvr$chrom <- gsub("chr", "", sapply(strsplit(cnvr$CNVR, ":"), '[', 1))
cnvr$Start <- as.numeric(sapply(strsplit(sapply(strsplit(cnvr$CNVR, ":"), "[", 2), "-"), "[", 1))
cnvr$End <- as.numeric(sapply(strsplit(sapply(strsplit(cnvr$CNVR, ":"), "[", 2), "-"), "[", 2))

### Read in SAIGE results for deletion-only and duplication-only model ================
dels <- read.table(paste(wkdir, "step2/DEL_SAIGE_step2_AN_rare.txt", sep="/"), header=TRUE)
dels$neglog10P <- -log10(dels$p.value)
dels$OR <- exp(dels$BETA)
dups <- read.table(paste(wkdir, "step2/DUP_SAIGE_step2_AN_rare.txt", sep="/"), header=TRUE)
dups$neglog10P <- -log10(dups$p.value)
dups$OR <- exp(dups$BETA)

### read in genes .bed file
genes <- read.table(paste(wkdir, "CNVR/Genes2", sep="/"))

## x-scale for plots
scaleFUN <- function(x) {
  tmp <- x/1000000
  return(tmp)
}

#=========== DUPLICATIONS =======================================================================
cnvrs <- cnvr %>% filter(Genotype=="Duplication")
cnvrs$N_Cases <- NA
cnvrs$N_Conts <- NA
cnvrs$N_Case_CNV <- NA
cnvrs$N_Cont_CNV <- NA

for (i in 1:nrow(cnvrs)){
  
  print(i)
  ## extract CNVR boundaries and a buffer region around the CNVR
  c <- cnvrs[i, "chrom"]
  start <- cnvrs[i, "Start"]
  end <- cnvrs[i, "End"]
  length <- (end-start)
  min <- start-length
  max <- end+length
  
  ## extract CNVs that intersect the CNVR
  cnvr_tmp <- c(c, start, end, "tmp")
  cnvr_tmp <- data.frame(t(unlist(cnvr_tmp)))
  fwrite(cnvr_tmp, "TMP1", col.names = F, sep = "\t")
  system(paste("plink --noweb --cfile /QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.DUP --cnv-intersect TMP1 --cnv-write --out TMP2"))
  tracks <- fread("TMP2.cnv")
  fam <- fread("TMP2.fam")
  tracks$status <- fam$V6[match(tracks$FID, fam$V1)]
  tracks$status <- ifelse(tracks$status==1, "Control", "Case")
  tracks$status <- as.factor(tracks$status)
  lower <- min(tracks$BP1)
  upper <- max(tracks$BP2)
  
  ## finding plot boundaries
  if (lower>min){
    xlim_lower <- min
  } else {
    xlim_lower <- lower
  }
  if (upper<max){
    xlim_upper <- max
  } else {
    xlim_upper <- upper
  }
  
  ### creating CNV track plot for cases
  cases <- tracks %>% filter(status=="Case")
  cnvrs$N_Case_CNV[i] <- nrow(cases)
  cnvrs$N_Cases[i] <- length(unique(cases$FID))
  cases$ID <- 1:nrow(cases)
  p_cases <- ggplot(cases, aes(x = BP1, y = ID)) +
    geom_segment(aes(xend = BP2, yend = ID), linewidth=3) +
    geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
    geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
    theme_classic(base_size = 100) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          legend.position = "none") + 
    labs(title="Case Duplications") +
    scale_x_continuous(labels =scaleFUN, limits=c(xlim_lower, xlim_upper))
  
  ### creating CNV track plot for controls
  controls <- tracks %>% filter(status=="Control")
  cnvrs$N_Cont_CNV[i] <- nrow(controls)
  cnvrs$N_Conts[i] <- length(unique(controls$FID))
  controls$ID <- 1:nrow(controls)
  p_controls <- ggplot(controls, aes(x = BP1, y = ID)) +
    geom_segment(aes(xend = BP2, yend = ID), linewidth=3) +
    geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
    geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
    theme_classic(base_size = 100) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none") + 
    labs(title="Control Duplications", x = paste0("Cytoband ", cnvrs[i, "Cytoband"], " (MB)")) +
    scale_x_continuous(labels =scaleFUN, limits=c(xlim_lower,xlim_upper))
  
  ## extract Saige results within CNVR with the corresponding genotype
  cnvr_saige <- dups %>% filter(CHR==c) %>% filter(MarkerID >= min & MarkerID <= max)
  cnvr_breakpoints <- unique(c(tracks$BP1, tracks$BP2))
  cnvr_saige <- cnvr_saige %>% filter(MarkerID%in%cnvr_breakpoints)
  
  ## breakpoint SAIGE plot
  cnvr_saige_steps <- cnvr_saige %>% filter(MarkerID >= start & MarkerID <= end)
  bplot <- ggplot() +
    geom_point(data=cnvr_saige, mapping=aes(x=MarkerID, y=neglog10P), size=15, color="black") +
    theme_classic(base_size=100) +
    geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
    geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
    geom_hline(yintercept=1.30103, linetype="dashed", linewidth=7, color="grey") +
    labs(y="-log10(P-value)",  
         title="SAIGE Breakpoint Associations (Duplication-only Model)") +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels =scaleFUN, limits=c(xlim_lower,xlim_upper))
  
  ## extract intersecting gene tracks and create plot
  if (!is.na(cnvrs[i, "Gene_Names"])){
    gene <- cnvrs[i, "Gene_Names"]
    gene <- strsplit(gene, " ")
    df_genes <- data.frame()
    for (g in 1:length(gene[[1]])){
      gene2 <- gene[[1]][g]
      bound <- genes[genes$V4==gene2,]
      bound$Y <- g
      df_genes <- rbind(df_genes, bound)
    }
    p_genes <- ggplot(df_genes, aes(x=V2, y=Y)) + 
      geom_segment(aes(x = V2, xend = V3, y = Y, yend = Y), linewidth = 5) +
      geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
      geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
      theme_classic(base_size = 100) +
      theme(legend.position="none",
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x=element_blank()) +
      labs(title="Gene Tracks") +
      ylim(0, 3) +
      scale_x_continuous(labels =scaleFUN, expand=c(0, 0), limits=c(xlim_lower,xlim_upper)) +
      geom_text_repel(aes(label = V4),
                      size = 20, color="black", max.overlaps=10)
  }
  
  ## all plots
  if (!is.na(cnvrs[i, "Gene_Names"])){
    png(filename = paste0(wkdir, "/CNVR/plots/", cnvrs[i, "chrom"], "_", start, "_",end, "_", cnvrs[i, "Genotype"], ".png"),
        width = 70, height = 60, units='in', bg="white", res=150, type=c("cairo"))
    print(bplot/p_genes/p_cases/p_controls +  plot_layout(heights = c(4,1,2,2)))
    dev.off()
  } else {
    png(filename = paste0(wkdir, "CNVR/plots/", cnvrs[i, "chrom"], "_", start, "_",end, "_", cnvrs[i, "Genotype"], ".png"),
        width = 70, height = 60, units='in', bg="white", res=150, type=c("cairo"))
    print(bplot/p_cases/p_controls +  plot_layout(heights = c(4,2,2)))
    dev.off()
  }
}



#==== DELETIONS =============================================================================================================
cnvrs <- cnvr %>% filter(Genotype=="Deletion")
cnvrs$N_Cases <- NA
cnvrs$N_Conts <- NA
cnvrs$N_Case_CNV <- NA
cnvrs$N_Cont_CNV <- NA


for (i in 1:nrow(cnvrs)){
  
  print(i)
  ## extract CNVR boundaries and some space outside of it
  c <- cnvrs[i, "chrom"]
  start <- cnvrs[i, "Start"]
  end <- cnvrs[i, "End"]
  length <- (end-start)
  min <- start-length
  max <- end+length
  
  ## extract CNVR boundaries and a buffer region around the CNVR
  cnvr_tmp <- c(c, start, end, "tmp")
  cnvr_tmp <- data.frame(t(unlist(cnvr_tmp)))
  fwrite(cnvr_tmp, "TMP1", col.names = F, sep = "\t")
  system(paste("plink --noweb --cfile /QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.DEL --cnv-intersect TMP1 --cnv-write --out TMP2"))
  tracks <- fread("TMP2.cnv")
  fam <- fread("TMP2.fam")
  tracks$status <- fam$V6[match(tracks$FID, fam$V1)]
  tracks$status <- ifelse(tracks$status==1, "Control", "Case")
  tracks$status <- as.factor(tracks$status)
  lower <- min(tracks$BP1)
  upper <- max(tracks$BP2)
  
  ## finding plot boundaries
  if (lower>min){
    xlim_lower <- min
  } else {
    xlim_lower <- lower
  }
  if (upper<max){
    xlim_upper <- max
  } else {
    xlim_upper <- upper
  }
  
  ### creating CNV track plot for cases
  cases <- tracks %>% filter(status=="Case")
  cnvrs$N_Case_CNV[i] <- nrow(cases)
  cnvrs$N_Cases[i] <- length(unique(cases$FID))
  cases$ID <- 1:nrow(cases)
  p_cases <- ggplot(cases, aes(x = BP1, y = ID)) +
    geom_segment(aes(xend = BP2, yend = ID), linewidth=3) +
    geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
    geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
    theme_classic(base_size = 100) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          legend.position = "none") + 
    labs(title="Case Deletions") +
    scale_x_continuous(labels =scaleFUN, limits=c(xlim_lower, xlim_upper))
  
  ### creating CNV track plot for controls
  controls <- tracks %>% filter(status=="Control")
  cnvrs$N_Cont_CNV[i] <- nrow(controls)
  cnvrs$N_Conts[i] <- length(unique(controls$FID))
  controls$ID <- 1:nrow(controls)
  p_controls <- ggplot(controls, aes(x = BP1, y = ID)) +
    geom_segment(aes(xend = BP2, yend = ID), linewidth=3) +
    geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
    geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
    theme_classic(base_size = 100) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none") + 
    labs(title="Control Deletions", x = paste0("Cytoband ", cnvrs[i, "Cytoband"], " (MB)")) +
    scale_x_continuous(labels =scaleFUN, limits=c(xlim_lower,xlim_upper))
  
  ## extract Saige results within CNVR with the corresponding genotype
  cnvr_saige <- dels %>% filter(CHR==c) %>% filter(MarkerID >= min & MarkerID <= max)
  cnvr_breakpoints <- unique(c(tracks$BP1, tracks$BP2))
  cnvr_saige <- cnvr_saige %>% filter(MarkerID%in%cnvr_breakpoints)
  
  ## breakpoint SAIGE plot
  cnvr_saige_steps <- cnvr_saige %>% filter(MarkerID >= start & MarkerID <= end)
  bplot <- ggplot() +
    geom_point(data=cnvr_saige, mapping=aes(x=MarkerID, y=neglog10P), size=15, colour="black") +
    theme_classic(base_size=100) +
    geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
    geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
    geom_hline(yintercept=1.30103, linetype="dashed", linewidth=7, color="grey") +
    labs(y="-log10(P-value)",  
         title="SAIGE Breakpoint Associations: (Deletion-only Model)") +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels =scaleFUN, limits=c(xlim_lower,xlim_upper))
  
  ## extract intersecting gene tracks and creating plot
  if (!is.na(cnvrs[i, "Gene_Names"])){
    gene <- cnvrs[i, "Gene_Names"]
    gene <- strsplit(gene, " ")
    df_genes <- data.frame()
    for (g in 1:length(gene[[1]])){
      gene2 <- gene[[1]][g]
      bound <- genes[genes$V4==gene2,]
      bound$Y <- g
      if (bound$V2 < xlim_lower){
        bound$V2 <- xlim_lower
      }
      if (bound$V3 > xlim_upper){
        bound$V3 <- xlim_upper
      }
      df_genes <- rbind(df_genes, bound)
    }
    
    p_genes <- ggplot(df_genes, aes(x=V2, y=Y)) + 
      geom_segment(aes(x = V2, xend = V3, y = Y, yend = Y), linewidth = 5) +
      geom_vline(xintercept=start, linetype="dashed", linewidth=7, color="grey") +
      geom_vline(xintercept=end, linetype="dashed", linewidth=7, color="grey") +
      theme_classic(base_size = 100) +
      theme(legend.position="none",
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x=element_blank()) +
      labs(title="Gene Tracks") +
      ylim(0, 3) +
      scale_x_continuous(labels =scaleFUN, expand=c(0,0), limits=c(xlim_lower,xlim_upper)) +
      geom_text_repel(aes(label = V4),
                      size = 20, color="black", max.overlaps=10)
  }
  
  ## all plots
  if (!is.na(cnvrs[i, "Gene_Names"])){
    png(filename = paste0(wkdir, "/CNVRS/plots/", cnvrs[i, "chrom"], "_", start, "_",end, "_", cnvrs[i, "Genotype"], ".png"),
        width = 70, height = 60, units='in', bg="white", res=150, type=c("cairo"))
    print(bplot/p_genes/p_cases/p_controls +  plot_layout(heights = c(4,1,2,2)))
    dev.off()
  } else {
    png(filename = paste0(wkdir, "CNVRS/plots/", cnvrs[i, "chrom"], "_", start, "_",end, "_", cnvrs[i, "Genotype"], ".png"),
        width = 70, height = 60, units='in', bg="white", res=150, type=c("cairo"))
    print(bplot/p_cases/p_controls +  plot_layout(heights = c(4,2,2)))
    dev.off()
  }
}









