geneMatrix_readme.txt

PF Sullivan, UNC/KI
updated 14 June 2022

Metadata for geneMatrix. Gene models from GENCODE v40 (4/2021), based on hg38
(but hg19 included). Created by 01geneMatrix.R. Tab-delimited text, strings
quoted, header plus 61498 rows, and 99 columns. Gene-based, all GENCODE GTF
genes. Coordinates are 1-based. 

order	column	pmid	variable_description
1	ensgid	34791404	ENSEMBL gene identifier stable (ENSG00000223972)
2	gene_name	34791404	gene name (STXBP1)
3	dup_gene_name	34791404	is gene_name duplicated?
4	gene_type	34791404	biotype: protein_coding, miRNA, lincRNA, etc
5	PAR	34791404	is gene in human pseudoautosomal region?
6	EMHC	34791404	is gene in extended MHC region? chr6:25-34 mb
7	gene_id	34791404	ENSEMBL gene identifier versioned (ENSG00000223972.11)
8	hg38h0	34791404	hg38 chromosome, chr1-22 chrX chrY chrM
9	h1	34791404	hg38 gene start
10	h2	34791404	hg38 gene end
11	hstr	34791404	hg38 gene strand
12	hg19g0	34791404	hg19 chromosome, chr1-22 chrX chrY chrM (may be missing or different)
13	g1	34791404	hg19 gene start
14	g2	34791404	hg19 gene end
15	gstr	34791404	hg19 gene strand
16	description	34791404	name of protein product
17	fracCdsPhylopAm	33177664	fraction of CDS constrained in 240 mammals, protein-coding, phyloP
18	fracCdsPhylopAm_decile	33177664	as above as decile (0-9)
19	fracCdsPhastConsPr	33177664	fraction of CDS constrained in 43 primates, protein-coding, phastCons
20	fracExonPhylopAm	33177664	fraction of exons (cDNA) constrained in 240 mammals, phyloP
21	fracExonPhastConsPr	33177664	fraction of exons (cDNA) constrained in 43 primates, phastCons
22	fracExonPhylopAm_decile	33177664	as above as decile (0-9)
23	loeuf	32461654	gnomAD, LOEUF, degree intolerance to pLoF (0=least)
24	loeuf_decile	32461654	gnomAD LOEUF decile (0=most constrained, 9=least)
25	pLI	32461654	gnomAD v2.1.1 pLI (0=least, 1=most)
26	NeQTLsnps	32913098	GTEx v8, number of unique eQTL-SNPs in any tissue
27	exprMeanDecile	32913098	GTEx v8, gene expression in body, decile (0=lowest, 9=highest)
28	GWAS	30445434	does this gene have any gwsig GWAS catalog hits ±10kb?
29	gwasCatalog	30445434	if TRUE, list of GWAS traits (separated by " // ")
30	OMIM	17357067	does this gene have an OMIM entry?
31	omimList	17357067	if TRUE, list of OMIM entries (separated by " // ")
32	mendelianPsychiatric	32122747	Mendelian with psychiatric symptoms?
33	dd2020.wes.fdr05	33057194	T/F Dev delay WES, Kaplanis Nature 2020
34	asd2020.wes.fdr05	31981491	T/F ASD WES, Satterstrom Cell 2020
35	asd2021.wes.fdr05	0	T/F ASD WES, Fu 2021 bioRxiv
36	dd2021.wes.fdr05	0	T/F Dev delay WES, Fu 2021 bioRxiv
37	ndd2021.wes.fdr05	0	T/F Neuro dev disorder WES, Fu 2021 bioRxiv
38	scz2022.wes.bonfSig05	35396579	T/F SCZ WES, SCHEMA, Nature 2022, bonferroni correction
39	scz2022.wes.qvalSig05	35396579	T/F SCZ WES, SCHEMA, Nature 2022, qvalue
40	IDgene	26748517	is this any type of ID gene? (sysID)
41	IDdominant	26748517	ID gene, dominant
42	IDrecessive	26748517	ID gene, recessive
43	IDhaploInsufficient	26748517	ID gene, haploinsufficient
44	CNVdelPheno	0	CNV deletion a=ASD d=dev delay i=ID s=SCZ
45	CNVdupPheno	0	CNV duplication a=ASD d=dev delay i=ID s=SCZ
46	synGOgene	31171447	SynGO gene, curated
47	inInversion	31269027	is this gene in a human inversion region?
48	tfgene	29425488	does this gene code a transcription factor?
49	housekeeping	23810203	is this a housekeeping gene? high expression, low CV over many tissues
50	entrezgene_id	0	Entrez gene ID
51	hgnc_id	0	HUGO gene nomenclature ID
52	asd2020.wes.qvalue	31981491	as above, but q-values not T/F
53	asd2021.wes.qvalue	0	as above, but q-values not T/F
54	dd2021.wes.qvalue	0	as above, but q-values not T/F
55	ndd2021.wes.qvalue	0	as above, but q-values not T/F
56	scz2022.wes.qvalue	35396579	as above, but q-values not T/F
57	scz2022InLocus	35396580	pgc scz2022 gwas, gene in LD-defined locus
58	scz2022PriorityGene	35396580	pgc scz2022 gwas, prioritized gene (fine mapping, eQTL, HiC, rare variant)
59	scz2022eqtl	35396580	pgc scz2022 gwas, eQTL-gene for eSNP in locus (GTEx cortex)
60	scz2022hic	35396580	pgc scz2022 gwas, Hi-C anchor in locus and anchor near gene TSS (HCRCI)
61	ensgid.mouse	0	ensembl gene id in mouse
62	mouseHuman1to1	0	ortholog, human-mouse, 1:1, high confidence
63	Adipose_Tissue	32913098	GTEx v8 mean expression TPM
64	Adrenal_Gland	32913098	GTEx v8 mean expression TPM
65	Blood	32913098	GTEx v8 mean expression TPM
66	Blood_Vessel	32913098	GTEx v8 mean expression TPM
67	BrainCortex	32913098	GTEx v8 mean expression TPM
68	BrainOther	32913098	GTEx v8 mean expression TPM
69	Breast	32913098	GTEx v8 mean expression TPM
70	Colon	32913098	GTEx v8 mean expression TPM
71	Esophagus	32913098	GTEx v8 mean expression TPM
72	Heart	32913098	GTEx v8 mean expression TPM
73	Liver	32913098	GTEx v8 mean expression TPM
74	Lung	32913098	GTEx v8 mean expression TPM
75	Muscle	32913098	GTEx v8 mean expression TPM
76	Nerve	32913098	GTEx v8 mean expression TPM
77	Ovary	32913098	GTEx v8 mean expression TPM
78	Pancreas	32913098	GTEx v8 mean expression TPM
79	Pituitary	32913098	GTEx v8 mean expression TPM
80	Prostate	32913098	GTEx v8 mean expression TPM
81	Salivary_Gland	32913098	GTEx v8 mean expression TPM
82	Skin	32913098	GTEx v8 mean expression TPM
83	Small_Intestine	32913098	GTEx v8 mean expression TPM
84	Spleen	32913098	GTEx v8 mean expression TPM
85	Stomach	32913098	GTEx v8 mean expression TPM
86	Testis	32913098	GTEx v8 mean expression TPM
87	Thyroid	32913098	GTEx v8 mean expression TPM
88	Uterus	32913098	GTEx v8 mean expression TPM
89	Vagina	32913098	GTEx v8 mean expression TPM
90	enstxid_reference	34791404	transcript selected for constrain and size measures (pickOne)
91	gene_size	0	gene boundary
92	transcript_size	0	transcript size
93	cdna_size	0	cDNA size, sum exon lengths (coding + noncoding)
94	cds_size	0	CDS size (protein-coding genes only)
95	SpeciesN	0	Number of species in ENSEMBL collection with 1:1 ortholog to human
96	SpeciesFrac	0	Fraction of species
97	PrimateN	0	Number of primates in ENSEMBL collection with 1:1 ortholog to human
98	PrimateFrac	0	Fraction of primates
99	exprMean	0	ignore


