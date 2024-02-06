
#!/bin/bash

# This script filters for rare CNVs (rCNVs) and splits by CNV type (dels and dups) and length before conducting any CNV burden analyses

#===first split hg38 .cnv file into chromosome-wide .cnv files====================================================

cfile="/QRISdata/Q4399/Anorexia/UKB/plink_files/UKBB_CNVs_for_AN_hg38"
output_dir="/QRISdata/Q4399/Anorexia/UKB/rare_cnvs"

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

output_dir="/QRISdata/Q4399/Anorexia/UKB/rare_cnvs"
path="/QRISdata/Q4399/Anorexia/UKB/plink_files"
file_name="UKBB_CNVs_for_AN_hg38"
i=$(awk 'BEGIN {FS=" "} END {print int((NR/100)+1)}' "${path}/${file_name}.fam")
p="/QRISdata/Q4399/software/plink-1.07-x86_64/plink"

echo -e "\tFID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES" > ${output_dir}/${file_name}.rare.cnv
touch ${output_dir}/${file_name}.rare.freq.cnv

for chr in {1..22}; do
  echo ${chr}
  cp ${path}/${file_name}.fam ${output_dir}/chr${chr}_${file_name}.fam
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${file_name} --cnv-make-map --out ${output_dir}/chr${chr}_${file_name}
  ## extract rare
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${file_name} --cnv-freq-exclude-above ${i} --cnv-write --out ${output_dir}/chr${chr}_${file_name}.rare
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${file_name}.rare --cnv-make-map --out ${output_dir}/chr${chr}_${file_name}.rare
  sed 1d ${output_dir}/chr${chr}_${file_name}.rare.cnv >> ${output_dir}/${file_name}.rare.cnv
  ## calculate cnv frequency
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${file_name}.rare --cnv-freq-method2 0.5 --cnv-write --cnv-write-freq --out ${output_dir}/chr${chr}_${file_name}.rare.freq
  sed 1d ${output_dir}/chr${chr}_${file_name}.rare.freq.cnv >> ${output_dir}/${file_name}.rare.freq.cnv
done

cp ${path}/${file_name}.fam ${output_dir}/${file_name}.rare.fam
${p} --noweb --cfile ${output_dir}/${file_name}.rare --cnv-make-map --out ${output_dir}/${file_name}.rare



#===================split rare CNVs into dup and del=============================================================

WKDIR="/QRISdata/Q4399/Anorexia/UKB/rare_cnvs"
file_name="UKBB_CNVs_for_AN_hg38"

## rare deletions
${p} --noweb --cfile ${WKDIR}/${file_name}.rare --cnv-del --cnv-write --out ${WKDIR}/${file_name}.DEL
${p} --noweb --cfile ${WKDIR}/${file_name}.DEL --cnv-make-map --out ${WKDIR}/${file_name}.DEL
## rare duplications
${p} --noweb --cfile ${WKDIR}/${file_name}.rare --cnv-dup --cnv-write --out ${WKDIR}/${file_name}.DUP
${p} --noweb --cfile ${WKDIR}/${file_name}.DUP --cnv-make-map --out ${WKDIR}/${file_name}.DUP
## all rare cnvs (deletions and duplications)
${p} --noweb --cfile ${WKDIR}/${file_name}.rare --cnv-write --out ${WKDIR}/${file_name}.ALL
${p} --noweb --cfile ${WKDIR}/${file_name}.ALL --cnv-make-map --out ${WKDIR}/${file_name}.ALL



#========= Split by length ======================================================================================

WKDIR="/QRISdata/Q4399/Anorexia/UKB/rare_cnvs"
file_name="UKBB_CNVs_for_AN_hg38"

for type in ALL DEL DUP; do

  echo ${type}
  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type} --cnv-kb 20 --cnv-max-kb 100 --cnv-write  --out ${WKDIR}/${file_name}.${type}.first
  ${p} --noweb --cfile ${WKDIR}/U${file_name}.${type}.first --cnv-make-map --out ${WKDIR}/${file_name}.${type}.first

  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type} --cnv-kb 100 --cnv-max-kb 200 --cnv-write  --out ${WKDIR}/${file_name}.${type}.second
  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type}.second --cnv-make-map --out ${WKDIR}/${file_name}.${type}.second

  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type} --cnv-kb 200 --cnv-max-kb 500 --cnv-write  --out ${WKDIR}/${file_name}.${type}.third
  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type}.third --cnv-make-map --out ${WKDIR}/${file_name}.${type}.third

  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type} --cnv-kb 500 --cnv-write  --out ${WKDIR}/${file_name}.${type}.fourth
  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type}.fourth --cnv-make-map --out ${WKDIR}/${file_name}.${type}.fourth

  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type} --cnv-write --out ${WKDIR}/${file_name}.${type}.all
  ${p} --noweb --cfile ${WKDIR}/${file_name}.${type}.all --cnv-make-map --out ${WKDIR}/${file_name}.${type}.all
done



