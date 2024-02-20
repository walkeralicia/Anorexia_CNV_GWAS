
#!/bin/bash
# 
#PBS -S /bin/bash
#PBS -A UQ-IMB
#PNS -N collins
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=10:mem=100GB

p="/QRISdata/Q4399/software/plink-1.07-x86_64/plink"
WKDIR="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs"
cnvdir="/QRISdata/Q4399/Anorexia/UKB/rare_cnvs"
file_name="UKBB_CNVs_for_AN_hg38"

# duplications ====================

for i in {1..101}
do
region=$(awk 'NR == '${i}' {print $4}' ${WKDIR}/collins_dups.txt)
echo ${region}
awk 'NR == '${i}' {print $1, $2, $3, $4}' ${WKDIR}/collins_dups.txt > ${WKDIR}/region.txt
mkdir -p ${WKDIR}/DUP/${region}
${p} --noweb --cfile ${cnvdir}/${file_name}.DUP --cnv-kb 100 --cnv-intersect ${WKDIR}/region.txt --cnv-overlap 0.5 --cnv-write --out ${WKDIR}/DUP/${region}/${file_name}.DUP.${region}
${p} --noweb --cfile ${WKDIR}/DUP/${region}/${file_name}.DUP.${region} --cnv-make-map --out ${WKDIR}/DUP/${region}/${file_name}.DUP.${region}

${p} --noweb --cfile ${cnvdir}/${file_name}.DUP --cnv-intersect ${WKDIR}/region.txt --cnv-region-overlap 0.25 --cnv-write --out ${WKDIR}/DUP/${region}/${file_name}.DUP.${region}.reciprocal
${p} --noweb --cfile ${WKDIR}/DUP/${region}/${file_name}.DUP.${region}.reciprocal --cnv-make-map --out ${WKDIR}/DUP/${region}/${file_name}.DUP.${region}.reciprocal
done


# deletions =================================

for i in {1..77}
do
region=$(awk 'NR == '${i}' {print $4}' ${WKDIR}/collins_dels.txt)
echo ${region}
awk 'NR == '${i}' {print $1, $2, $3, $4}' ${WKDIR}/collins_dels.txt > ${WKDIR}/region.txt
mkdir -p ${WKDIR}/DEL/${region}
${p} --noweb --cfile ${cnvdir}/${file_name}.DEL --cnv-kb 100 --cnv-intersect ${WKDIR}/region.txt --cnv-overlap 0.5 --cnv-write --out ${WKDIR}/DEL/${region}/${file_name}.DEL.${region}
${p} --noweb --cfile ${WKDIR}/DEL/${region}/${file_name}.DEL.${region} --cnv-make-map --out ${WKDIR}/DEL/${region}/${file_name}.DEL.${region}

${p} --noweb --cfile ${cnvdir}/${file_name}.DEL --cnv-intersect ${WKDIR}/region.txt --cnv-region-overlap 0.25 --cnv-write --out ${WKDIR}/DEL/${region}/${file_name}.DEL.${region}.reciprocal
${p} --noweb --cfile ${WKDIR}/DEL/${region}/${file_name}.DEL.${region}.reciprocal --cnv-make-map --out ${WKDIR}/DEL/${region}/${file_name}.${region}.reciprocal
done

