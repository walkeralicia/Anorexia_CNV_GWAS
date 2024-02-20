# Genome-wide CNV Association Study

 This Github contains a pipeline to conduct a genome wide copy number variant association study (CNV-GWAS) and accompanies the article "Genome-wide CNV-association study in 8,674 individuals with anorexia nervosa". 

 # Citations

If you use scripts from this Github please cite the article: ref (TBC). 
 
 # Pipeline

Folders **/01_Samples** and **/02_CNVs** describe Anorexia Nervosa (AN) phenotype extraction and CNV-level quality-control specific to acquired UKB input data. 

All downstream folders describe a general CNV-burden pipeline that can be applied on any CNV data that is in PLINK cfile format (i.e., .cnv, .fam, .map). Moreover, any set of covariates (e.g., PCs, Sex, Age, Array type, etc.) and trait (binary or continuous) of interest can be used. 

 Folders with scripts:

-  **/01_Samples** extracts the UKB AN phenotypes.

- **/02_CNVs** processes the acquired UKB CNV calls into PLINK cfile format.

- **/03_CNV_Burden** extracts rare (<1% population frequency) CNVs and contains subfolders:

     **/genome_wide** to conduct various total genome-wide CNV burden analyses.
   
     **/locus_wide** to conduct multiple locus-wide CNV burden analyses using two sets of CNV lists; a set of 167 dosage-sensitive pleiotropic CNVs; and a set of 67 well-established syndromic CNVs.

- **/04_CNV_breakpoint_GWAS** converts CNV breakpoint input data from PLINK cfile format into PLINK bfile format to then conduct a CNV-GWAS.

- **/05_Novel_CNV_regions** identifies novel disease-risk CNV regions (CNVRs) and plots the results.

- **/06_Meta_analyses** meta-analyses locus-wide and CNV-GWAS results with a replication study's results using Stouffer's method. The subfolder **/ANGI_data** contains ANGI summary statistics used in the study to meta-analyse with the UKB.

# CNV Input Data

 In this pipeline, CNV input data is acquired from the UK Biobank (UKB) called by Kendal et al. Due to privacy restrictions, UKB CNV-GWAS results are meta-analysed with summary statistics from the Anorexia Nervosa Genetics Initiaitve (ANGI) study using Stouffer's method. 

# CNV Annotation Input Data

CNV annotation files come from various public databases and published research papers. All files and corresponding README files have been placed in the folder **/data**. When using this data, please cite the sources appropriately. 

1. **Syndromic CNV List**
2. **Dosage-sensitive CNV List**
3. **Zoonomia scores**
4. **Gene annotation Files**


 

 

