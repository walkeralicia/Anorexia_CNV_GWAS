

#===This script generates  UKB CNV-breakpoint frequency plots (not included within published paper)=======

#=== load required R packages
library(data.table)
library(dplyr)
library(ggplot2)
library(egg)
library(ggrepel)

#========= Read in SAIGE results for CNV breakpoints of any frequency============================
input <- "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE/step2"
output <- "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE/plots"
dels <- read.table(paste(input, "DEL_SAIGE_step2_AN.txt", sep="/"), header=T)
dups <- read.table(paste(input, "DUP_SAIGE_step2_AN.txt", sep="/"), header=T)

#======= Read in Genomic cytoband information===============================================
cytoband <- fread(paste(output, "cytoBand.txt", sep="/"), header=FALSE) %>% as.data.frame()
cytoband$CHR <- as.numeric(gsub("chr", "", cytoband$V1))

#======= calculating cnv frequency using genotype counts ====================================
fam <- read.table("/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.fam")
numpeople <- nrow(fam)
dels2 <- dels %>% select(CHR, MarkerID, AC_Allele2, AF_Allele2)
dups2 <- dups %>% select(CHR, MarkerID, AC_Allele2, AF_Allele2)
both <- merge(dels2, dups2, by = c("CHR", "MarkerID"), all=TRUE)
both$NumCNV <- both$AC_Allele2.y + both$AC_Allele2.x
both$DUP_CNVF <- (100*both$AC_Allele2.y)/numpeople
both$DEL_CNVF <- (100*both$AC_Allele2.x)/numpeople
both$CNVF <- (100*both$NumCNV)/numpeople

#============ ordering probes ===============================================================
both$Group <- ifelse(both$CHR%in%c(1,3,5,7,9,11,13,15,17,19,21), "odd", "even")
both_ord <- both[order(both$CHR, both$MarkerID, decreasing=FALSE),]
both_ord$ID <- 1:nrow(both_ord)


#====== mapping probes to cytobands=========================================================
both_ord$cytoband <- NA
for (i in 1:dim(both_ord)[1]){
  print(i)
  chr = both_ord[i, "CHR"]
  pos = both_ord[i, "MarkerID"]
  cytoband_chr = cytoband[cytoband$CHR==chr,]
  cytoband_chr$pos = ifelse(cytoband_chr$V2 <= pos & cytoband_chr$V3 >= pos, "yes", "no")
  cyto <- cytoband_chr[which(cytoband_chr$pos=="yes"), "V4"]
  both_ord$cytoband[i] <- cyto
}

#========Finding top probes for each cytoband ==========================================
both_ord$CHR_cytoband <- paste0("chr", both_ord$CHR, both_ord$cytoband)
uni_cyto <- unique(both_ord$CHR_cytoband)
both_ord$topdupprobe <- NA
both_ord$topdelprobe <- NA
for (i in 1:length(uni_cyto)){
  print(uni_cyto[i])
  both_ord_cyto <- both_ord[both_ord$CHR_cytoband==uni_cyto[i],]
  chr=both_ord_cyto$CHR[1]
  top_dup_probe <- both_ord_cyto[which(both_ord_cyto$DUP_CNVF==max(both_ord_cyto$DUP_CNVF, na.rm=TRUE)), "MarkerID"][1]
  dup_pos <- which(both_ord$MarkerID==top_dup_probe & both_ord$CHR==chr)
  both_ord$topdupprobe[dup_pos] <- "yes"
  top_del_probe <- both_ord_cyto[which(both_ord_cyto$DEL_CNVF==max(both_ord_cyto$DEL_CNVF, na.rm=TRUE)), "MarkerID"][1]
  del_pos <- which(both_ord$MarkerID==top_del_probe & both_ord$CHR==chr)
  both_ord$topdelprobe[del_pos] <- "yes"
}


#===== Highlight and label probes that have a cnv frequency of greater than 3 =======================
both_ord$highlightdup <- ifelse(!is.na(both_ord$topdupprobe) & both_ord$DUP_CNVF > 3, "yes", "no")
both_ord$highlightdel <- ifelse(!is.na(both_ord$topdelprobe) & both_ord$DEL_CNVF > 3, "yes", "no")
both_ord$labeldup <- ifelse(both_ord$highlightdup=="yes", both_ord$CHR_cytoband, NA)
both_ord$labeldel <- ifelse(both_ord$highlightdel=="yes", both_ord$CHR_cytoband, NA)

#==== Save results so far ==========================================================================================
write.table(both_ord, paste(output,"cnv_freq_plot.txt", sep="/"), col.names=T, row.names=F, sep="\t", quote=F)


#=====generate CNV-b frequency plot for duplications==========================================
p1 <- ggplot(data=both_ord, aes(x=ID, y=DUP_CNVF, color=Group)) +
  geom_point(size=4)  +
  geom_text_repel(aes(label = labeldup),
                  size = 10, color="black") + 
  labs(x='', y = 'Duplications') +
  scale_color_manual(values = c("odd" = "#0099FF", "even" = "333333")) + 
  theme_classic(base_size=50) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none") +
  ggtitle("UKB CNV Frequency (%)") +
  xlim(0, nrow(both_ord))

png(filename =  paste(output,"dup_freq.png", sep="/"),width = 30, height = 15, units='in', bg="white", res=150, type=c("cairo"))
print(p1)
dev.off()

#======generate CNV-b frequency plot for deletions=======================================================                            
p2 <- ggplot(data=both_ord, aes(x=ID, y=DEL_CNVF, color=Group)) +
  geom_point(size=4)  +
  geom_text_repel(aes(label = labeldel),
                  size = 10, color="black") + 
  labs(x='', y = 'Deletions') +
  scale_color_manual(values = c("odd" = "#FF9933", "even" = "#333333")) + 
  theme_classic(base_size=50) +
  scale_y_reverse() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none") +
  xlim(0, nrow(both_ord))

png(filename =  paste(output,"del_freq.png", sep="/"),width = 30, height = 15, units='in', bg="white", res=150, type=c("cairo"))
print(p2)
dev.off()    

#=========== combine the plots ============================================================================
plots <- ggarrange(p1, ggplot() + theme_void(),
                   p2 + 
                     theme(axis.title.x = element_blank() ), 
                   nrow=3, 
                   heights = c(1, -0.1, 1))

png(filename =  paste(output,"cnv_freq.png", sep="/"),width = 34, height = 15, units='in', bg="white", res=150, type=c("cairo"))
print(plots)
dev.off()
                            