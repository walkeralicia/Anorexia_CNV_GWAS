# Genome-wide CNV Association Study

 This Github contains a pipeline to conduct a genome-wide copy-number-variant association study (CNV-GWAS) and accompanies the article "Genome-wide copy number variant association study with anorexia nervosa". 

 # Citations

If you use scripts from this Github please consider citing the article: ref (TBC). 
 
 # Pipeline

All folders downstream of and including **/03_CNV_Burden** describe a CNV burden pipeline that can be applied on any CNV data that is in PLINK cfile format (i.e., .cnv, .fam, .map). Moreover, the pipeline can be applied to any binary or continuous trait of interest and to any set of covariates (e.g., principal components, sex, age, array type).

Folders **/01_Samples** and **/02_CNVs** describe Anorexia Nervosa (AN) phenotype extraction and CNV-level quality-control specific to acquired UK Biobank (UKB) input data. 

 Folders with scripts:

-  **/01_Samples** extracts the UKB AN phenotypes.

- **/02_CNVs** processes the acquired UKB CNV calls into PLINK cfile format.

- **/03_CNV_Burden** extracts rare (<1% population frequency) CNVs and contains subfolders:

     **/S2_genome_wide** to conduct total genome-wide CNV burden analyses.
   
     **/S2_locus_wide** to conduct locus-wide CNV burden analyses using two sets of CNV lists; 1) a set of 178 dosage-sensitive pleiotropic CNVs; 2) a set of 67 well-established syndromic CNVs.

- **/04_CNV_breakpoint_GWAS** converts CNV breakpoint input data from PLINK cfile format into PLINK bfile format to then conduct a CNV-GWAS.

- **/05_Novel_CNV_regions** identifies novel disease-risk CNV regions (CNVRs) and plots results.

- **/06_Meta_analyses** meta-analyses locus-wide and CNV-GWAS results with a replication study's results using Stouffer's method. The subfolder **/ANGI_data** contains Anorexia Nervosa Genetics Initiative (ANGI) summary statistics used in the cited study to meta-analyse with the UKB.

# CNV Input Data

In this pipeline, CNV input data is acquired from the UKB called by Kendal et al. Due to privacy restrictions, UKB CNV-GWAS results are meta-analysed with summary statistics from ANGI using Stouffer's method. 

# CNV Annotation Input Data

CNV annotation files come from various public databases and published research papers. All files and corresponding README files have been placed in the folder **/data**. When using this data, please cite the sources appropriately. 

1. **Syndromic CNV List**
2. **Dosage-sensitive CNV List**
3. **Zoonomia scores**
4. **Gene annotation Files**


 

 

