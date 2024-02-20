
#!/bin/bash

#= This script filters for rare CNVs (rCNVs) and splits by CNV type (dels and dups) and length ===================================

#== Set up paths for script=======================================================================

in_dir="/Anorexia/UKB/plink_files"
output_dir="/Anorexia/UKB/rare_cnvs"
cfile="UKBB_CNVs_for_AN_hg38"

#=== Split hg38 .cnv file into chromosome-wide .cnv files ========================================

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through chromosomes and extract "rare" CNVs
for i in {1..22}; do
    echo "Processing Chromosome $i"
    
    # Extract CNVs for the current chromosome
    awk -v chr="$i" '$3 == chr' "$in_dir/${cfile}.cnv" > "$output_dir/chr${i}_${cfile}.cnv"
    
    # Print the number of CNVs for the current chromosome
    echo "Number of CNVs for Chromosome $i: $(wc -l < "$output_dir/chr${i}_${cfile}.cnv")"
done



#=== Apply Frequency(<1%) filter for each Chr and combine genome-wide================================

i=$(awk 'BEGIN {FS=" "} END {print int((NR/100)+1)}' "${in_dir}/${cfile}.fam")
p="plink"

echo -e "\tFID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES" > ${output_dir}/${cfile}.rare.cnv
touch ${output_dir}/${cfile}.rare.freq.cnv

for chr in {1..22}; do
  
  echo ${chr}
  cp ${in_dir}/${cfile}.fam ${output_dir}/chr${chr}_${cfile}.fam
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${cfile} --cnv-make-map --out ${output_dir}/chr${chr}_${cfile}
  
  ## extract rare
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${cfile} --cnv-freq-exclude-above ${i} --cnv-write --out ${output_dir}/chr${chr}_${cfile}.rare
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${cfile}.rare --cnv-make-map --out ${output_dir}/chr${chr}_${cfile}.rare
  sed 1d ${output_dir}/chr${chr}_${cfile}.rare.cnv >> ${output_dir}/${cfile}.rare.cnv
  
  ## calculate cnv frequency
  ${p} --noweb --cfile ${output_dir}/chr${chr}_${cfile}.rare --cnv-freq-method2 0.5 --cnv-write --cnv-write-freq --out ${output_dir}/chr${chr}_${cfile}.rare.freq
  sed 1d ${output_dir}/chr${chr}_${cfile}.rare.freq.cnv >> ${output_dir}/${cfile}.rare.freq.cnv
done

cp ${in_dir}/${cfile}.fam ${output_dir}/${cfile}.rare.fam
${p} --noweb --cfile ${output_dir}/${cfile}.rare --cnv-make-map --out ${output_dir}/${cfile}.rare



#===================split rare CNVs into dup and del=============================================================

## rare deletions
${p} --noweb --cfile ${output_dir}/${cfile}.rare --cnv-del --cnv-write --out ${output_dir}/${cfile}.DEL
${p} --noweb --cfile ${output_dir}/${cfile}.DEL --cnv-make-map --out ${output_dir}/${cfile}.DEL

## rare duplications
${p} --noweb --cfile ${output_dir}/${cfile}.rare --cnv-dup --cnv-write --out ${output_dir}/${cfile}.DUP
${p} --noweb --cfile ${output_dir}/${cfile}.DUP --cnv-make-map --out ${output_dir}/${cfile}.DUP

## all rare cnvs (deletions and duplications)
${p} --noweb --cfile ${output_dir}/${cfile}.rare --cnv-write --out ${output_dir}/${cfile}.ALL
${p} --noweb --cfile ${output_dir}/${cfile}.ALL --cnv-make-map --out ${output_dir}/${cfile}.ALL



#========= Split by length ======================================================================================


for type in ALL DEL DUP; do

  echo ${type}
  # 20-100kb
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type} --cnv-kb 20 --cnv-max-kb 100 --cnv-write  --out ${output_dir}/${cfile}.${type}.first
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type}.first --cnv-make-map --out ${output_dir}/${cfile}.${type}.first
  # 100-200kb
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type} --cnv-kb 100 --cnv-max-kb 200 --cnv-write  --out ${output_dir}/${cfile}.${type}.second
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type}.second --cnv-make-map --out ${output_dir}/${cfile}.${type}.second
  # 200-500kb
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type} --cnv-kb 200 --cnv-max-kb 500 --cnv-write  --out ${output_dir}/${cfile}.${type}.third
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type}.third --cnv-make-map --out ${output_dir}/${cfile}.${type}.third
  #>500kb
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type} --cnv-kb 500 --cnv-write  --out ${output_dir}/${cfile}.${type}.fourth
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type}.fourth --cnv-make-map --out ${output_dir}/${cfile}.${type}.fourth
  # all
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type} --cnv-write --out ${output_dir}/${cfile}.${type}.all
  ${p} --noweb --cfile ${output_dir}/${cfile}.${type}.all --cnv-make-map --out ${output_dir}/${cfile}.${type}.all

done



