# CNV Genome-wide Association Study

 This Github contains a workflow pipeline to conduct a copy number variant genome wide association study (CNV-GWAS) and accompanies the article "Genome-wide CNV-association study in 8,674 individuals with anorexia nervosa". 

 # Citations

If you use scripts from this Github please cite the article: ref (TBC). 

 # CNV Input Data

 In this pipeline, CNV input data is acquired from the UK Biobank (UKB) called by Kendal et al. Due to privacy restrictions, UKB CNV-GWAS results are meta-analysed with summary statistics from the Anorexia Nervosa Genetics Initiaitve (ANGI) study using Stouffer's method. 
 
 # Pipeline

Folders **/01_Samples** and **/02_CNVs** describe AN phenotype extraction and CNV-level quality-control specific to acquired UKB input data. 

All downstream folders describe a general CNV-burden pipeline that can be applied on any CNV data that is in PLINK cfile format (i.e., .cnv, .fam, .map) and on any set of covariates (e.g., PCs, Sex, Age, Array type, etc.). 

 Folders with scripts:

-  **/01_Samples** extracts the UKB Anorexia Nervosa (AN) phenotypes.

- **/02_CNVs** processes the acquired UKB CNV calls into PLINK cfile format.

- **/03_CNV_Burden** extracts rare (<1% population frequency) CNVs and contains subfolders:

     **/genome_wide** to conduct various total genome-wide CNV burden analyses.
   
     **/locus_wide** to conduct multiple locus-wide CNV associations using two sets of CNV lists; a set of 167 dosage-sensitive, pleiotropic CNVs, and a set of 67 well-established syndromic CNVs.

- **/04_CNV_breakpoint_GWAS** converts CNV input data from PLINK cfile format into PLINK bfile format to then conduct a CNV-GWAS.

- **/05_Novel_CNV_regions** identifies novel disease-risk CNV regions (CNVRs) and plots the results.

- **/06_Meta_analyses_with_ANGI** meta-analyses locus-wide and CNV-GWAS results with a replication study's results using Stouffer's method. The subfolder **/ANGI_data** contains ANGI summary statistics used in the study to meta-analyse with the UKB.

# CNV Annotation Input Data

CNV annotation files come from various public databases and published research papers. All files and corresponding README files have been placed in the folder **/data**. When using this data, please cite the sources appropriately. 

1. **Syndromic CNV List**
2. **Dosage-sensitive CNV List**
3. **Zoonomia scores**
4. **Gene annotation Files**


 

 

