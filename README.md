# Genome-wide CNV Association Study

 This Github contains the workflow pipeline to conduct a rare Copy Number Variant (rCNV) genome-wide association study (GWAS) within the UK Biobank (UKB) using Anorexia Nervosa (AN) and Body Mass Index (BMI) as the outcomes. The pipeline also meta-analyses UKB summary statistics with the Anorexia Nervosa Genetics Initiative (ANGI) study using Stouffer's method. 

 # Citations

 If you use scripts from this Github please consider citing the article: 

 # Pipeline

 Folders **/01_Samples** and **/02_CNVs** describe AN phenotype extraction and CNV-level quality-control specific to acquired UKB input data. All downstream folders describe a general CNV-burden workflow pipeline that can be applied on any CNV data that is in PLINK cfile format (i.e., .cnv, .fam, .map).

 Folders with scripts:

-  **/01_Samples** extracts the UKB Anorexia Nervosa (AN) phenotypes.

- **/02_CNVs** processes the acquired UKB CNV calls into PLINK cfile format.

- **/03_CNV_Burden** extracts rare (<1% population frequency) CNVs (rCNVs) and contains subfolders:

     **/genome_wide** to conduct various total genome-wide rCNV burden analyses.
   
     **/locus_wide** to conduct multiple locus-wide rCNV associations using two sets of CNV lists; a set of 167 dosage-sensitive, pleiotropic CNVs, and a set of 67 well-established syndromic CNVs.

- **/04_CNV_breakpoint_GWAS** converts cnv input data from PLINK cfile format into PLINK bfile format to then conduct a rCNV breakpoint genome-wide association study (rCNV-GWAS).

- **/05_Novel_CNV_regions** identifies novel disease-risk CNV regions (CNVRs) and plots the results.

- **/06_Meta_analyses_with_ANGI** meta-analyses locus-wide and rCNV-breakpoint GWAS association results with a replication study's results using Stouffer's method. The subfolder **/ANGI_data** contains ANGI summary statistics used in the study to meta-analyse with the UKB.

 # CNV Input Data

 In this study, UKB CNV input data comes from the acquired CNVs called on autosomes by Kendal et al.  All UKB samples were genotyped using two Affymetrix arrays (UK BiLEVE and UK Biobank Axiom arrays), each containing over 800,000 probes. Kendal et al. used PennCNV to identify all CNVs that spanned at least 10 probes and were longer than 20 kb. 

# CNV Annotation Input Data

CNV annotation files come from various public databases and published research papers. All files and corresponding README files have been placed in the folder **/data**. When using this data, please cite the sources appropriately. 

1. **Syndromic CNV List**
2. **Dosage-sensitive CNV List**
3. **Zoonomia scores**
4. **Gene annotation Files**


 

 

