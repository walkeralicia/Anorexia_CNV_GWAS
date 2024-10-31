# Genome-wide CNV Association Study

 This Github contains a pipeline to conduct a genome-wide copy-number-variant association study (CNV-GWAS) and accompanies the article "Genome-wide copy number variation association study in anorexia nervosa" published in Molecular Psychiatry. 

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

- **/06_Meta_analyses** meta-analyses locus-wide results with a replication study's results using Stouffer's method. The subfolder **/ANGI_data** contains Anorexia Nervosa Genetics Initiative (ANGI) summary statistics for the 21 novel AN-associated CNVR used in the cited study to meta-analyse with the UKB.

# CNV Input Data

In this pipeline, CNV input data was acquired from the UKB called by Kendal et al (DOI: 10.1192/bjp.2018.301). Due to privacy restrictions, UKB CNV-GWAS results were meta-analysed with summary statistics from the Anorexia Nervosa Genetics Initiative (ANGI) using Stouffer's method. ANGI is a multi-country collaboration that includes research teams in the US, Sweden, Denmark, and Australia with assistance from New Zealand. ANGI's objective is to identify the genetic causes of AN. Previous publications (Thornton et al., DOI: 10.1016/j.cct.2018.09.015, and Watson et al., DOI: 10.1038/s41588-019-0439-2) have provided comprehensive information on participant recruitment, phenotyping, DNA collection and genotyping. CNVs were identified using Illumina Global Screening Array raw intensity data (618,540 probes) from 13,787 individuals. Autosomal and X chromosomal CNVs were identified using EnsemblemCNV, a software wrapper that incorporates three CNV calling algorithms: iPattern, PennCNV, and QuantiSNP. CNVs were considered valid if they were detected by a minimum of 2 algorithms, had a length of more than 20kb, included at least 10 probes, and had successfully passed a series of quality control measures outlined in the supplementary materials of the publication associated with this Github ("Genome-wide copy number variation association study in anorexia nervosa").

# CNV Annotation Input Data

CNV annotation files come from various public databases and published research papers. All files have been placed in the folder **/data**. When using this data, please cite the sources appropriately. 

1. **Syndromic CNV List**

We compiled a list of 67 syndromic developmental-associated CNVs (**CNV_List.xlsx**), including 45 deletions and 22 duplications, to test for association with AN status. Of these 67 CNVs, 60 were extracted from the DECIPHER database v11.19 (DOI: 10.1016/j.ajhg.2009.03.010). Selection of the remaining 7 syndromic CNVs were detailed in the methods of the publication associated with this Github ("Genome-wide copy number variation association study in anorexia nervosa").
 
2. **Dosage-sensitive CNV List**

We tested a CNV list of 178 dosage-sensitive genomic segments (**1-s2.0-S0092867422007887-mmc3.xlsx** and **1-s2.0-S0092867422007887-mmc4.xlsx**), encompassing 77 deletions and 101 duplications, that confer disease-risk across 54 complex and Mendelian traits/disorders, including 24 neurodevelopmental traits. We further assessed genome-wide CNV burden in AN as a function of CNVs intersecting haploinsufficient or triplo-sensitive protein-coding genes, using scores provided within **collins2022-TableS7-pHaplo-pTriplo.tsv**. These three lists were sourced from Collins et al. (DOI: 10.1016/j.cell.2022.06.036). 

3. **Zoonomia scores**

Evolutionary constraint scores, obtained from the Zoonomia Consortium (DOI: 10.1038/s41586-020-2876-6), were used to functionally annotate all novel AN-associated CNV regions, and to assess genome-wide CNV burden with AN as a function of the average proportion of bases within CNVs that are evolutionary constrained. We could not provide the raw zoonomia scores, although we have provided the associated README file (**zoonomia-readme.docx**). 
  
4. **Other annotation files**

The files **cytoBand.txt** and **NRXN1_exons.csv** provide genomic structural information and were both sourced from the UCSC Genome Browser. The file **geneMatrix.tsv** provides information on genes from GENCODE v40 (4/2021) with metadata found within **geneMatrix_readme.txt**. These files were kindly provided by Patrick Sullivan.
 

 

