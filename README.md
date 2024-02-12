# Genome-wide CNV Association Study

 This Github contains the workflow pipeline to conduct a Genome Wide Copy Number Variant (CNV) association study within the UK Biobank (UKB) using Anorexia Nervosa (AN) as the outcome. 

 # Citations

 If you use scripts from this Github please consider citing the article: 

 # Pipeline

 Folders 01 and 02 are UKB specific, describing AN phenotype extraction and CNV-level QC. However, all subsequent folders describe a general CNV-burden analyses pipeline that can be applied on any CNV input data that is in PLINK cfile format (i.e., .cnv, .fam, .map). 

1. Folder /01_Samples contains the script to extract UKB Anorexia Nervosa (AN) phenotypes.
2. Folder 02/CNVs contains the script to process the acquired UKB CNV calls into PLINK cfile format (i.e., .cnv, .map, .fam files) for further CNV burden analyses.
3. Folder /03_CNV_Burden contains the R script (S1_extract_rare_CNVs.R) to extract rare (<1% population frequency) CNVs and contains:
   a. Subfolder /01_genome_wide: contains scripts to conduct total rare genome-wide CNV burden anlayses, split by CNV-type, CNV length, CNV count, CNV frequency, proportion of mammalian constraint bases covered, and number of haploinsufficient or triplosensitive genes.
   b. Subfolder /02_locus_wide: contains scripts to conduct locus-wide CNV associations studies using two sets of known CNV lists; a set of 167 pleiotropic dosage-sensitive disease-risk CNVs, and a set of XX DECIPHER genomic syndrome CNVs.
4. Folder /04_CNV_breakpoint_GWAS: 


 # CNV Input Data

 In this Study, CNV input data comes from the acquired CNVs called on autosomes within the UKB by Kendal et al.  All UKB samples were genotyped using two Affymetrix arrays (UK BiLEVE and UK Biobank Axiom arrays), each containing over 800,000 probes. Kendal et al. used PennCNV to identify all CNVs that spanned at least 10 probes and were longer than 20 kb. 
 
# CNV Annotation Input data

CNV annotation files come from various public databases and published research papers, and have been placed in the folder /data. When using this data, please cite the sources appropriately. 

1. Developmental CNV List
2. Pleiotropic Disease-risk dosage-sensitive CNV List
3. Zoonomia scores
4. Gene annotations


 

 

