# Genome-wide CNV Association Study

 This Github contains the workflow pipeline to conduct a Genome Wide Copy Number Variant (CNV) association study within the UK Biobank (UKB) using Anorexia Nervosa (AN) and Body Mass Index (BMI) as the outcomes.

 # Citations

 If you use scripts from this Github please consider citing the article: 

 # Pipeline

 Folders **/01_Samples** and **/02_CNVs** describe AN phenotype extraction and CNV-level QC specific to acquired input data from the UKB. All subsequent folders describe a general CNV-burden workflow pipeline that can be applied on any CNV input data that is in the required PLINK cfile format (i.e., .cnv, .fam, .map).

 Folders with scripts:

-  **/01_Samples** extracts the UKB Anorexia Nervosa (AN) phenotypes.

- **/02_CNVs** processes the acquired UKB CNV calls into PLINK cfile format (i.e., .cnv, .map, .fam files) for further CNV burden analyses.

- **/03_CNV_Burden** extracts rare (<1% population frequency) CNVs and contains subfolders:

   a. **/01_genome_wide** to conduct various total genome-wide rCNV burden analyses.
   
   b. **/02_locus_wide** to conduct multiple locus-wide CNV associations using two sets of known CNV lists; a set of 167 pleiotropic dosage-sensitive disease-risk CNVs, and a set of XX DECIPHER genomic syndrome CNVs.

- **/04_CNV_breakpoint_GWAS** converts cnv input data from PLINK cfile format into PLINK bfile format. Downstream scripts ...

- **/05_novel_CNV_regions** identifies novel disease-associated CNV regions (CNVRs) and plots these CNVRs with genomic annotations.

- **/06_Meta_analyses_with_ANGI** meta-analyses locus-wide and CNV-breakpoint GWAS association results with a replication study's results using Stouffer's method. The folder also contains the summary statistics from the ANGI study used to meta-analyse with the UKB.

 # CNV Input Data

 In this Study, CNV input data comes from the acquired CNVs called on autosomes within the UKB by Kendal et al.  All UKB samples were genotyped using two Affymetrix arrays (UK BiLEVE and UK Biobank Axiom arrays), each containing over 800,000 probes. Kendal et al. used PennCNV to identify all CNVs that spanned at least 10 probes and were longer than 20 kb. 
 
# CNV Annotation Input data

CNV annotation files come from various public databases and published research papers, and have been placed in the folder /data. When using this data, please cite the sources appropriately. 

1. Developmental CNV List
2. Pleiotropic Disease-risk dosage-sensitive CNV List
3. Zoonomia scores
4. Gene annotations


 

 

