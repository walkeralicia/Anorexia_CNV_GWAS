
#========= This script creates plink bfiles from cfiles=================================================

#========== R script (cnv_to_plink_chrCHROMOSOME_batchBATCH.R)=======================

### READ IN DATA  
path<-"/QRISdata/Q4399/Anorexia/UKB/plink_files" # set your path to plink cfiles here
cnvs <- read.table(paste(path, "UKBB_CNVs_for_AN_hg38.cnv", sep="/"), header=TRUE)
fam <- read.table(paste(path, "UKBB_CNVs_for_AN_hg38.fam", sep="/"), header=FALSE)
num<-signif(nrow(fam),3)/1000

## Read in required r packages
library(dplyr)
library(data.table)
library(gtools)

## subset by chromosome
cnvs_chr <- cnvs[cnvs$CHR==CHROMOSOME, ]
breakpoints <- unique(c(cnvs_chr$BP1, cnvs_chr$BP2))
samples=fam$V1
batch = BATCH
subset_U <- round(length(samples)/num)*batch
subset_L <- subset_U - 999
if (subset_U == num*1000){
  subset_U = length(samples)
}
c(subset_L, subset_U)
samples = samples[subset_L:subset_U]

## map CNVs to corresponding allelic encoding
all_samples_dup <- list()
all_samples_del <- list()
for (i in 1:length(samples)){
  print(i)
  matrix_dup = matrix(data=NA, nrow=1, ncol=length(breakpoints))
  matrix_del = matrix(data=NA, nrow=1, ncol=length(breakpoints))
  cnvs_chr_sample <- cnvs_chr[cnvs_chr$FID == samples[i], ]
  if (is.na(cnvs_chr_sample[1,1])){
    vec <- rep("A A", length(breakpoints))
    matrix_dup <- matrix(vec, ncol = length(breakpoints), byrow = FALSE)
    matrix_del <- matrix(vec, ncol = length(breakpoints), byrow = FALSE)
  } else {
    for(j in 1:length(breakpoints)){
      cnvs_chr_sample_bp <- cnvs_chr_sample %>%
        mutate(cnv = breakpoints[j] >= cnvs_chr_sample$BP1 & breakpoints[j] <= cnvs_chr_sample$BP2)
      cnvs_chr_sample_bp <- cnvs_chr_sample_bp[cnvs_chr_sample_bp$cnv==TRUE,]
      cnvs_chr_sample_bp <- cnvs_chr_sample_bp[1,]
      if (!samples[i]%in%cnvs_chr_sample_bp$FID){
        matrix_dup[1,j] <- "A A"
        matrix_del[1,j] <- "A A"
      } else if (cnvs_chr_sample_bp[cnvs_chr_sample_bp["FID"]==samples[i], "TYPE"] < 2){
        matrix_dup[1,j] <- "0 0"
        matrix_del[1,j] <- "A T"
      } else {
        matrix_dup[1,j] <- "A T"
        matrix_del[1,j] <- "0 0"
      }
    }
  }
  all_samples_dup[[i]] <- matrix_dup
  all_samples_del[[i]] <- matrix_del
}

final_dup <- do.call(rbind, all_samples_dup)
colnames(final_dup) <- breakpoints
rownames(final_dup) <- samples
final_del <- do.call(rbind, all_samples_del)
colnames(final_del) <- breakpoints
rownames(final_del) <- samples


## .fam files
fam_batch <- fam[fam$V1%in%samples,]
colnames(fam_batch)<- c("FID", "IID","PID", "MID", "Sex", "Pheno")
rownames(fam_batch) <- samples
fam_batch <- fam_batch[match(rownames(final_dup), fam_batch$FID),]
fwrite(fam_batch, paste(path, "cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH.fam", sep="/"), col.names = F, row.names = F, quote = F, sep = "\t")
fwrite(fam_batch, paste(path, "cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH.fam", sep="/"), col.names = F, row.names = F, quote = F, sep = "\t")

## .map files
plink_ped_batch_dup <- cbind(fam_batch, as.data.frame(final_dup))
probe_order_dup <- colnames(plink_ped_batch_dup)[7:ncol(plink_ped_batch_dup)]
plink_map_dup <- data.frame(CHR = CHROMOSOME, SNP = probe_order_dup, GD = 0, POS = 1:length(probe_order_dup))
fwrite(plink_map_dup, paste(path, "cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH.map", sep="/"), col.names = F, row.names = F, quote = F, sep = "\t")
plink_ped_batch_del <- cbind(fam_batch, as.data.frame(final_del))
probe_order_del <- colnames(plink_ped_batch_del)[7:ncol(plink_ped_batch_del)]
plink_map_del <- data.frame(CHR = CHROMOSOME, SNP = probe_order_del, GD = 0, POS = 1:length(probe_order_del))
fwrite(plink_map_del, paste(path, "cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH.map", sep="/"), col.names = F, row.names = F, quote = F, sep = "\t")


# Order variants according to the .map file, split each variant column into two columns (to mimic a .ped file), and save the file.
plink_ped_batch_v2_dup <- plink_ped_batch_dup[, c(colnames(plink_ped_batch_dup)[1:6], probe_order_dup)]
plink_ped_batch_v2_dup <- cbind(plink_ped_batch_v2_dup[,1:6], as.data.frame(unlist(lapply(plink_ped_batch_v2_dup[, c(7:ncol(plink_ped_batch_v2_dup))], data.table::tstrsplit, " "), recursive = FALSE)))
fwrite(plink_ped_batch_v2_dup, paste(path, "cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH.ped", sep="/"), col.names = F, row.names = F, quote = F, sep = "\t")
plink_ped_batch_v2_del <- plink_ped_batch_del[, c(colnames(plink_ped_batch_del)[1:6], probe_order_del)]
plink_ped_batch_v2_del <- cbind(plink_ped_batch_v2_del[,1:6], as.data.frame(unlist(lapply(plink_ped_batch_v2_del[, c(7:ncol(plink_ped_batch_v2_del))], data.table::tstrsplit, " "), recursive = FALSE)))
fwrite(plink_ped_batch_v2_del, paste(path, "cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH.ped", sep="/"), col.names = F, row.names = F, quote = F, sep = "\t")

system("path=/QRISdata/Q4399/Anorexia/UKB/plink_files")
system("mkdir -p ${path}/bfiles")
system("plink --noweb --file ${path}/cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH --make-bed --out ${path}/bfiles/cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH")
system("rm ${path}/cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH.map")
system("rm ${path}/cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH.ped")
system("rm ${path}/cnv_to_PLINK_DUP_chrCHROMOSOME_batchBATCH.fam")
system("plink --noweb --file ${path}/cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH --make-bed --out ${path}/bfiles/cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH")
system("rm ${path}/cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH.map")
system("rm ${path}/plink_files/cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH.ped")
system("rm ${path}/plink_files/cnv_to_PLINK_DEL_chrCHROMOSOME_batchBATCH.fam")

