
#!/bin/bash


# This script runs S1_cnv_to_plink_chrCHROMOSOME_batchBATCH.R

#============= generate job script to run  S1_cnv_to_plink_chrCHROMOSOME_batchBATCH.R ============================

echo '#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N chrCHROMOSOME_cnv_to_plink
#PBS -A UQ-IMB
#PBS -l walltime=24:00:00
#PBS -l select=1:mem=50GB
#PBS -J 1-387

module load R/4.1.2
export R_LIBS=/QRISdata/Q4399/R_libraries/new_rlib_4.1.2

path="/QRISdata/Q4399/Anorexia/UKB/plink_files"

R --file=${path}/jobs/S1_cnv_to_plink_chrCHROMOSOME_batch${PBS_ARRAY_INDEX}.R' > ${path}/S1_cnv_to_plink_chrCHROMOSOME.sh

#================= Submit job script======================================================================

path="/QRISdata/Q4399/Anorexia/UKB/plink_files"
num=$(awk 'BEGIN {FS=" "} END {print int(NR/1000)}' "${path}/UKBB_CNVs_for_AN_hg38.fam")

for i in {1..22}; do
  sed 's/CHROMOSOME/'${i}'/g' ${path}/S1_cnv_to_plink_chrCHROMOSOME_batchBATCH.R > ${path}/S1_cnv_to_plink_chr${i}_batchBATCH.R
  sed 's/CHROMOSOME/'${i}'/g' ${path}/S1_cnv_to_plink_chrCHROMOSOME.sh > ${path}/S1_cnv_to_plink_chr${i}.sh

  j=1
  while [ $j -le $num ]; do
    sed 's/BATCH/'${j}'/g' ${path}/S1_cnv_to_plink_chr${i}_batchBATCH.R > ${path}/S1_cnv_to_plink_chr${i}_batch${j}.R
    ((j++))
  done

  qsub ${path}/S1_cnv_to_plink_chr${i}.sh
done


