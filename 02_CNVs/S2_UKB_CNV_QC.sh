#!/bin/bash

# This script creates plink .map cfile, performs hg19 to hg38 genomic co-ordinate liftover, 
# and generates an additional set of cfiles with BMI as the UKB phenotypic outcome

#============ create plink .map file for cfile data  (PLINK V1.07)===================================
# File paths
p="plink"
cfile="/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN"

${p} --noweb --cfile ${cfile} --cnv-make-map --out ${cfile}


#=============== hg19 to hg38 liftover ===============================================================
# liftover software and chain file
liftOver="liftOver"
chain_file="hg19ToHg38.over.chain"

# Perform .cnv to .bed conversion
# .bed format is chr, BP1, BP2, ID
awk 'BEGIN {OFS="\t"} NR>1 {print "chr" $3, $4, $5, NR-1}' "${cfile}.cnv" > "${cfile}_preliftover.bed"

# Perform liftover
${liftOver} ${cfile}_preliftover.bed ${chain_file} ${cfile}_postliftover.bed ${cfile}_postliftover_unMapped


#================= Create new set of plink cfiles with hg38 coordinates=========================================
# First add ID column to .cnv file
awk 'BEGIN {OFS="\t"} NR>1 {print NR-1, $0}' "${cfile}.cnv" > "${cfile}_temp.cnv"

# Join files based on the ID column
IFS=$'\t' # Set IFS to tab
echo -e "FID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES" > "${cfile}_hg38.cnv"
# Sort the input files
sort -k1,1 "${cfile}_temp.cnv" > "${cfile}_temp.sorted.cnv"
sort -k4,4 "${cfile}_postliftover.bed" > "${cfile}_postliftover.sorted.bed"
join -1 1 -2 4 -o 1.2,1.3,1.4,2.2,2.3,1.7,1.8,1.9 "${cfile}_temp.sorted.cnv" "${cfile}_postliftover.sorted.bed" | tr ' ' '\t' >> "${cfile}_hg38.cnv"
unset IFS 

# Clean up temporary file
rm "${cfile}_temp.cnv" "${cfile}_temp.sorted.cnv" "${cfile}_postliftover.sorted.bed" 


#================ create new plink .fam and .map cfiles with hg38 coordinates=========================================
cp ${cfile}.fam ${cfile}_hg38.fam
${p} --noweb --cfile ${cfile}_hg38 --cnv-make-map --out ${cfile}_hg38


#=============== create a second set of plink cfiles for BMI as the outcome ============================================
# File paths
cfile="/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38"
bmi_file="OP_year_080623_participant.csv"
map_file="ukb12505bridge14421.txt"
output_file='${cfile}"_BMI"'

# Step 1: Extract relevant columns from the fam file
awk '{print $2, $3, $4, $5, $6}' "${cfile}.fam" > "${cfile}.temp.fam"

# Step 2: Extract relevant columns from the map file
awk '{print $1, $2}' "${map_file}" > "${output_file}.temp.map"

# Step 3: Extract relevant columns from the BMI file
awk -F',' 'NR>1{print $1, $2}' "${bmi_file}" > "${output_file}.temp.bmi"

# Step 4: Join the temp.map and temp.bmi files based on the first column
awk 'NR==FNR{a[$1]=$2; next} {print $1, $2, a[$1]}' "${output_file}.temp.bmi" "${output_file}.temp.map" > "${output_file}.temp.join1"

# Step 5: Join the temp.join1 and temp.fam file
awk 'NR==FNR{a[$1]=$3; next} {print $1, $1, $2, $3, $4, a[$1]}' "${output_file}.temp.join1" "${cfile}.temp.fam" > "${output_file}.temp.fam.join1"

# Clean up temporary files
rm "${cfile}.temp.fam" "${output_file}.temp.map" "${output_file}.temp.bmi" "${output_file}.temp.join1"
mv "${output_file}.temp.fam.join1" ${output_file}_temp.fam

# replace empty cells for BMI with -9
awk -v OFS="\t" '{if ($6=="") $6=-9; print}' ${output_file}_temp.fam > ${output_file}.fam

# create plink cfiles (.cnv and .map for BMI outcome)
cp ${cfile}.cnv ${output_file}.cnv
${p} --noweb --cfile ${output_file} --cnv-make-map --out ${output_file}



