# Genome-wide CNV Association Study

 This Github contains a pipeline to conduct a genome-wide copy-number-variant association study (CNV-GWAS) and accompanies the article "Genome-wide copy number variation association study in anorexia nervosa". 

 # Citations

If you use scripts from this Github please consider citing the article: ref (TBC). 
 
 # Pipeline

All folders from **/03_CNV_Burden**  to **/06_Meta_analyses** describe a CNV burden pipeline that can be applied on any CNV data that is in PLINK cfile format (i.e., .cnv, .fam, .map).

Folders **/01_Samples** and **/02_CNVs** describe Anorexia Nervosa (AN) phenotype extraction and CNV-level quality-control specific to acquired UK Biobank (UKB) input data. 

 Folders with scripts:

-  **/01_Samples** extracts the UKB AN phenotypes.

- **/02_CNVs** processes the acquired UKB CNV calls into PLINK cfile format.

- **/03_CNV_Burden** extracts rare (<1% population frequency) CNVs and contains subfolders:

     **/S2_genome_wide** to conduct total genome-wide CNV burden analyses.
   
     **/S2_locus_wide** to conduct locus-wide CNV burden analyses using two sets of CNV lists; 1) a set of 178 dosage-sensitive pleiotropic CNVs; 2) a set of 67 well-established syndromic CNVs.

- **/04_CNV_breakpoint_GWAS** converts CNV breakpoint input data from PLINK cfile format into PLINK bfile format to then conduct a CNV-GWAS.

- **/05_Novel_CNV_regions** identifies novel disease-risk CNV regions (CNVRs) and plots results.

- **/06_Meta_analyses** meta-analyses locus-wide results with a replication study's results using Stouffer's method. The subfolder **/ANGI_data** contains Anorexia Nervosa Genetics Initiative (ANGI) summary statistics used in the cited study to meta-analyse with the UKB. This encompasses summary statistics for CNVs associated with AN status, BMI in AN cases, and BMI in AN controls. 

# CNV Input Data

In this pipeline, CNV input data was acquired from the UKB called by Kendal et al. Due to privacy restrictions, UKB CNV-GWAS results were meta-analysed with summary statistics from the Anorexia Nervosa Genetics Initiative (ANGI) using Stouffer's method. ANGI is a multi-country collaboration that includes research teams in the US, Sweden, Denmark, and Australia with assistance from New Zealand. ANGI's objective is to identify the genetic causes of AN. Previous publications (Thornton et al. and Watson et al.) have provied comprehensive information on participant recruitment, phenotyping, DNA collection and genotyping. CNVs were identified using Illumina Global Screening Array raw intensity data (618,540 probes) from 13,787 individuals. Autosomal and X chromosomal CNVs were identified using EnsemblemCNV, a software wrapper that incorporates three CNV calling algorithms: iPattern, PennCNV, and QuantiSNP. CNVs were considered valid if they were detected by a minimum of 2 algorithms, had a length of more than 20kb, included at least 10 probes, and had successfully passed a series of quality control measures outlined in the supplementary materials of the associated publication ("Genome-wide copy number variation association study in anorexia nervosa").

# CNV Annotation Input Data

CNV annotation files come from various public databases and published research papers. All files have been placed in the folder **/data**. When using this data, please cite the sources appropriately. 

1. **Syndromic CNV List**
   - 
3. **Dosage-sensitive CNV List**
4. **Zoonomia scores**
5. **Gene annotation Files**


 

 

