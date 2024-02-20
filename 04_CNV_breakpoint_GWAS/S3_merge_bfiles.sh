
#!/bin/bash

#======= This script combines all chromosome-wide bfiles into a genome-wide bfile split by CNV type (deletions and duplications) =============

p="plink"
path="/Anorexia/UKB/plink_files"
outdir="${path}/bfiles"

for i in {1..22}; do
  echo ${i}
  touch ${outdir}/cnv_to_PLINK_DEL_chr${i}.list
  touch ${outdir}/cnv_to_PLINK_DUP_chr${i}.list
  j=1
  while [ $j -le $num ]; do
    echo -e cnv_to_PLINK_DEL_chr${i}_batch${j} >> ${outdir}/cnv_to_PLINK_DEL_chr${i}.list
    echo -e cnv_to_PLINK_DUP_chr${i}_batch${j} >> ${outdir}/cnv_to_PLINK_DUP_chr${i}.list
    ((j++))
  done
  ${p} --noweb --merge-list ${outdir}/cnv_to_PLINK_DEL_chr${i}.list --make-bed --out ${outdir}/cnv_to_PLINK_DEL_chr${i}
  ${p} --noweb --merge-list ${outdir}/cnv_to_PLINK_DUP_chr${i}.list --make-bed --out ${outdir}/cnv_to_PLINK_DUP_chr${i}
done


#================= merge all chromosomes===========================================================================
touch ${outdir}/cnv_to_PLINK_DEL.list
touch ${outdir}/cnv_to_PLINK_DUP.list

for c in {1..22}; do
  echo -e cnv_to_PLINK_DEL_chr${c} >> ${outdir}/cnv_to_PLINK_DEL.list
  echo -e cnv_to_PLINK_DUP_chr${c} >> ${outdir}/cnv_to_PLINK_DUP.list
done

${p} --noweb --merge-list ${outdir}/cnv_to_PLINK_DEL.list --make-bed --out ${outdir}/cnv_to_PLINK_DEL
${p} --noweb --merge-list ${outdir}/cnv_to_PLINK_DUP.list --make-bed --out ${outdir}/cnv_to_PLINK_DUP
