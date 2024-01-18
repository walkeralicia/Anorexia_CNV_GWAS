
#!/bin/bash

#========= This script runs CNV-breakpoint GWAS for all and only rare CNV breakpoints=============================================


#===First create .pheno file for runnning SAIGE-GWAS==================================

# Set paths
pheno_file="/QRISdata/Q4399/Anorexia/UKB/pheno/All_pheno_for_UKBB_filtered.dat"
fam_file="/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38.fam"
bmi_file="/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38_BMI.fam"
output_path="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE"
output_file="${output_path}/UKBB_CNVs_for_AN_hg38_SAIGE.pheno"

# R script content
script=$(cat <<EOF
pheno <- read.table("$pheno_file", header=TRUE)
fam <- read.table("$fam_file")
BMI <- read.table("$bmi_file")
fam\$V5 <- factor(fam\$V5)
fam\$V6 <- factor(fam\$V6)
fam\$Array <- pheno\$Array[match(fam\$V1, pheno\$V1)]
fam\$Sex <- pheno\$Sex[match(fam\$V1, pheno\$V1)]
fam\$Age <- pheno\$Age[match(fam\$V1, pheno\$V1)]
fam\$AN <- ifelse(fam\$V6 == 2, 1, 0)
fam\$BMI <- BMI\$V6[match(fam\$V1, BMI\$V1)]

write.table(fam, "$output_file", col.names=T, row.names=F, sep="\t", quote=F)
EOF
)


# Save R script to a file
echo "$script" > ${output_path}/generate_pheno.R

# Run R script using Rscript
Rscript ${output_path}/generate_pheno.R


#================Step 1 SAIGE for duplications (repeat for DELS)====================================

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --job-name=DUP_SAIGE_step1_AN
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o ${outdir}/step1/DUP_SAIGE_step1_AN.stdout
#SBATCH -e ${outdir}/DUP_SAIGE_step1_AN.stderr

export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

SAIGE=/home/uqawal15/SAIGE/extdata
bfile=/QRISdata/Q4399/Anorexia/UKB/plink_files/bfiles
outdir=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE/

Rscript ${SAIGE}/step1_fitNULLGLMM.R \
--plinkFile=${bfile}/cnv_to_PLINK_DUP \
--phenoFile=${outdir}/UKBB_CNVs_for_AN_hg38_SAIGE.pheno \
--phenoCol=AN \
--covarColList=V5,Array,Age \
--qCovarColList=V5,Array \
--sampleIDColinphenoFile=V1 \
--traitType=binary \
--outputPrefix=${outdir}/step1/DUP_SAIGE_step1_AN \
--nThreads=24 \
--IsOverwriteVarianceRatioFile=TRUE' > ${outdir}/DUP_SAIGE_step1_AN.sh

sbatch ${outdir}/DUP_SAIGE_step1_AN.sh



#============STEP 2 SAIGE for duplications (repeat for deletions)============================

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --job-name=DUP_SAIGE_step2_AN
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o ${outdir}/step2/DUP_SAIGE_step2_AN.stdout
#SBATCH -e ${outdir}/step2/DUP_SAIGE_step2_AN.stderr

export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

SAIGE=/home/uqawal15/SAIGE/extdata
bfile=/QRISdata/Q4399/Anorexia/UKB/plink_files/bfiles
outdir=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE/

Rscript ${SAIGE}/step2_SPAtests.R \
--bedFile=${bfile}/cnv_to_PLINK_DUP.bed \
--bimFile=${bfile}/cnv_to_PLINK_DUP.bim \
--famFile=${bfile}/cnv_to_PLINK_DUP.fam \
--AlleleOrder=alt-first \
--SAIGEOutputFile=${outdir}/step2/DUP_SAIGE_step2_AN.txt \
--minMAF=0 \
--minMAC=1 \
--GMMATmodelFile=${outdir}/step1/DUP_SAIGE_step1_AN.rda \
--varianceRatioFile=${outdir}/step1/DUP_SAIGE_step1_AN.varianceRatio.txt \
--is_Firth_beta=TRUE \
--pCutoffforFirth=0.05 \
--is_output_moreDetails=TRUE \
--LOCO=FALSE' > ${outdir}/DUP_SAIGE_step2_AN.sh

sbatch ${outdir}/DUP_SAIGE_step2_AN.sh


#==========rare variants========================================================================

bfile="/QRISdata/Q4399/Anorexia/UKB/plink_files/bfiles"
p="/QRISdata/Q4399/software/plink"

### filter for rare variants (plink v1.9) 
${p}  --bfile ${bfile}/cnv_to_PLINK_DEL --max-maf 0.01 --make-bed --out ${bfile}/cnv_to_PLINK_DEL.rare
${p}  --bfile ${bfile}/cnv_to_PLINK_DUP --max-maf 0.01 --make-bed --out ${bfile}/cnv_to_PLINK_DUP.rare

########## SAIGE step 2 for rare variants (duplications - repeat for deletions)########################

echo '#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --job-name=DUP_SAIGE_step2_AN_rare
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o ${outdir}/step2/DUP_SAIGE_step2_AN_rare.stdout
#SBATCH -e ${outdir}/step2/DUP_SAIGE_step2_AN_rare.stderr

export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

SAIGE=/home/uqawal15/SAIGE/extdata
bfile=/QRISdata/Q4399/Anorexia/UKB/plink_files/bfiles
outdir=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/SAIGE/

Rscript ${SAIGE}/step2_SPAtests.R \
--bedFile=${bfile}/cnv_to_PLINK_DUP.update.rare.bed \
--bimFile=${bfile}/cnv_to_PLINK_DUP.update.rare.bim \
--famFile=${bfile}/cnv_to_PLINK_DUP.update.rare.fam \
--AlleleOrder=alt-first \
--SAIGEOutputFile=${outdir}/step2/DUP_SAIGE_step2_AN_rare.txt \
--minMAF=0 \
--minMAC=1 \
--GMMATmodelFile=${outdir}/step1/DUP_SAIGE_step1_AN.rda \
--varianceRatioFile=${outdir}/step1/DUP_SAIGE_step1_AN.varianceRatio.txt \
--is_Firth_beta=TRUE \
--pCutoffforFirth=0.05 \
--is_output_moreDetails=TRUE \
--LOCO=FALSE' > ${outdir}/DUP_SAIGE_step2_AN_rare.sh

sbatch ${outdir}/DUP_SAIGE_step2_AN_rare.sh

