# Genome-wide CNV Association Study

 This Github contains the workflow pipeline to conduct a Genome Wide Copy Number Variant (CNV) association study within the UK Biobank (UKB) using Anorexia Nervosa (AN) and Body Mass Index (BMI) as the outcomes.

 # Citations

 If you use scripts from this Github please consider citing the article: 

 # Pipeline

 Folders **/01_Samples** and **/02_CNVs** describe AN phenotype extraction and CNV-level QC specific to acquired input data from the UKB. All subsequent folders describe a general CNV-burden workflow pipeline that can be applied on any CNV input data that is in the required PLINK cfile format (i.e., .cnv, .fam, .map).

-  **/01_Samples** contains the script to extract UKB Anorexia Nervosa (AN) phenotypes.

- **/02_CNVs** contains the script to process the acquired UKB CNV calls into PLINK cfile format (i.e., .cnv, .map, .fam files) for further CNV burden analyses.

- **/03_CNV_Burden** contains the script (S1_extract_rare_CNVs.R) to extract rare (<1% population frequency) CNVs (rCNV) and contains subfolders:

   a. **/01_genome_wide** with scripts to conduct total genome-wide rCNV burden analyses, split by CNV type (i.e., duplication or deletion), CNV length (ranging from 20kb to >500kb), CNV count, and CNV frequency. The scripts also conduct CNV-burden analyses investigating the proportion of mammalian constraint bases covered, and the number of intersecting dosage-sensitive, haploinsufficient or triplosensitive, genes.
   
   b. **/02_locus_wide** with scripts to conduct multiple locus-wide CNV associations studies using two sets of known CNV lists; a set of 167 pleiotropic dosage-sensitive disease-risk CNVs, and a set of XX DECIPHER genomic syndrome CNVs.

- **/04_CNV_breakpoint_GWAS** contains scripts to convert PLINK cfile format into PLINK bfile format with all CNV breakpoints encoded as a variant. Downstream scripts ...

- **/05_novel_CNV_regions** contains scripts to identify novel disease-associated CNV regions (CNVRs) and to generate plots of these CNVRs annotated with CNV-breakpoint GWAS association results, intersecting gene tracks, and intersecting CNV case-control tracks.

- **/06_Meta_analyses_with_ANGI** contains scripts to meta-analyse locus-wide and CNV-breakpoint GWAS association results with a replication study's results (e.g. ANGI) using Stouffer's method. The folder also contains the summary statistics from the ANGI study. 

 # CNV Input Data

 In this Study, CNV input data comes from the acquired CNVs called on autosomes within the UKB by Kendal et al.  All UKB samples were genotyped using two Affymetrix arrays (UK BiLEVE and UK Biobank Axiom arrays), each containing over 800,000 probes. Kendal et al. used PennCNV to identify all CNVs that spanned at least 10 probes and were longer than 20 kb. 
 
# CNV Annotation Input data

CNV annotation files come from various public databases and published research papers, and have been placed in the folder /data. When using this data, please cite the sources appropriately. 

1. Developmental CNV List
2. Pleiotropic Disease-risk dosage-sensitive CNV List
3. Zoonomia scores
4. Gene annotations


 

 

