# Genome-wide CNV Association Study

 This Github contains the workflow pipeline to conduct a genome-wide Copy Number Variant (CNV) association study within the UK Biobank (UKB) using Anorexia Nervosa (AN) as the outcome. 

 # Citations

 If you use scripts from this Github please consider citing the article: 

 # Pipeline

1. Folder /01_Samples contains the script to extract UKB Anorexia Nervosa (AN) phenotypes.
2. Folder 02/CNV_QC contains the script to process the UKB CNV calls into PLINK cfile format (i.e., .cnv, .map, .fam files) for further CNV burden analyses.

In this Study, CNV input data comes from the UKB, which can only be accessed through application. However, any CNV input once formatted into the required PLINK cfile format can be used in the downstream pipeline. 

3. Folder /03_CNV_burden contains the scripts to extract rare CNVs (<1% population frequnecy) and to perform lo



# CNV Annotation Input data

CNV annotation files come from various public databases and published research papers, and have been placed in the folder /data. When using this data, please cite the sources appropriately. 

1. Developmental CNV List
2. Pleiotropic Disease-risk dosage-sensitive CNV List
3. Zoonomia scores
4. Gene annotations


 

 

