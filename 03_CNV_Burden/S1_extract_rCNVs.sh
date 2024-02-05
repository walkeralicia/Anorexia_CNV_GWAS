
#!/bin/bash

# This script filters for rare CNVs and splits by CNV type (dels and dups) and length before conducting any CNV burden analyses

#===first split hg38 .cnv file into chromosome-wide .cnv files====================================================

cfile="/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38"
output_dir="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden/rare/"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through chromosomes and extract "rare" CNVs
for i in {1..22}; do
    echo "Processing Chromosome $i"
    
    # Extract CNVs for the current chromosome
    awk -v chr="$i" '$3 == chr' "$cfile.cnv" > "$output_dir/chr${i}_UKBB_CNVs_for_AN_hg38.cnv"
    
    # Print the number of CNVs for the current chromosome
    echo "Number of CNVs for Chromosome $i: $(wc -l < "$output_dir/chr${i}_UKBB_CNVs_for_AN_hg38.cnv")"
done



#===APPLY FREQUENCY (<1%) filter for each Chr and combine genome-wide========================================

output_dir="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden/rare/"
path="/QRISdata/Q4399/Anorexia/UKB/plink_files"
i=$(awk 'BEGIN {FS=" "} END {print int((NR/100)+1)}' "${path}/UKBB_CNVs_for_AN_hg38.fam")
p="/QRISdata/Q4399/software/plink-1.07-x86_64/plink"

echo -e "FID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES" > ${output_dir}/UKBB_CNVs_for_AN_hg38.rare.cnv
touch ${output_dir}/UKBB_CNVs_for_AN_hg38.rare.freq.cnv

for chr in {1..22}; do
  echo ${chr}
  cp ${path}/UKBB_CNVs_for_AN_hg38.fam ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.fam
  ${p} --noweb --cfile ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38 --cnv-make-map --out ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38
  ${p} --noweb --cfile ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38 --cnv-freq-exclude-above ${i} --cnv-write --out ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare
  ${p} --noweb --cfile ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare --cnv-make-map --out ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare
  sed 1d ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare.cnv >> ${output_dir}/UKBB_CNVs_for_AN_hg38.rare.cnv
  ${p} --noweb --cfile ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare --cnv-freq-method2 0.5 --cnv-write --cnv-write-freq --out ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare.freq
  sed 1d ${output_dir}/chr${chr}_UKBB_CNVs_for_AN_hg38.rare.freq.cnv >> ${output_dir}/UKBB_CNVs_for_AN_hg38.rare.freq.cnv
done

cp ${cfile}.fam ${output_dir}/UKBB_CNVs_for_AN_hg38.rare.fam
${p} --noweb --cfile ${output_dir}/UKBB_CNVs_for_AN_hg38.rare --cnv-make-map --out ${output_dir}/UKBB_CNVs_for_AN_hg38.rare

#===================split rare CNVs into dup and del=============================================================

WKDIR="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden"

${p} --noweb --cfile ${WKDIR}/rare/UKBB_CNVs_for_AN_hg38.rare --cnv-del --cnv-write --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.DEL
${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.DEL --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.DEL
${p} --noweb --cfile ${WKDIR}/rare/UKBB_CNVs_for_AN_hg38.rare --cnv-dup --cnv-write --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.DUP
${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.DUP --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.DUP
cp ${WKDIR}/rare/UKBB_CNVs_for_AN_hg38.rare.cnv ${WKDIR}/UKBB_CNVs_for_AN_hg38.ALL.cnv
cp ${WKDIR}/rare/UKBB_CNVs_for_AN_hg38.rare.fam ${WKDIR}/UKBB_CNVs_for_AN_hg38.ALL.fam
${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.ALL --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.ALL



#========= Split by length ======================================================================================

for type in ALL DEL DUP; do

  echo ${type}
  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type} --cnv-kb 20 --cnv-max-kb 100 --cnv-write  --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.first
  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.first --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.first

  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type} --cnv-kb 100 --cnv-max-kb 200 --cnv-write  --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.second
  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.second --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.second

  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type} --cnv-kb 200 --cnv-max-kb 500 --cnv-write  --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.third
  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.third --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.third

  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type} --cnv-kb 500 --cnv-write  --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.fourth
  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_filter.${type}.fourth --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.fourth

  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type} --cnv-write --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.all
  ${p} --noweb --cfile ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.all --cnv-make-map --out ${WKDIR}/UKBB_CNVs_for_AN_hg38.${type}.all
done



