# Genotype-variants
An R script for the visualization of RNA-seq variants from genetically engineered mouse models

What is Genotype variants.R?

Genotype Variants is a R script for adequate visualization of RNA-seq variants on the mouse mm10 build, based on VCF files (Variant Calling Format files, http://samtools.github.io/hts-specs/VCFv4.3.pdf). This R script takes as input two filtered VCF files (from wild-type and KO/KI genotypes, respectively) and output genome-wide plots of variants per genotype, including a summary of chromosomes with differential distribution of variants among genotypes. 

Input Preparation

To collect RNA-seq variants, BAM files from the alignment of RNA-seq reads are sorted and genome-wide simple diploid calling is performed by using Freebayes (https://github.com/ekg/freebayes). Raw VCF files are further processed using the VCFlib toolkit (https://github.com/vcflib/vcflib). VCF files needs to be intersected among biological replicates from each genotype (VCF-VCF intersect tool) and the resulting VCF file is then filtered using and adequate criteria. Finally, multiple allelic primitives (gaps or mismatches) are splitted using the Vcfallelicprimitives tool from the VCFlib toolkit. This pipeline can be fully implemented using the galaxy tools (https://usegalaxy.org/) setting the adequate version of the tools, as follows (May 11, 2018):

1) Input: Sorted BAM file from alignment (e.g HISAT)
2) FreeBayes: Simple diploid calling, default settings (Galaxy Version 1.1.0.46-0)
3) VCFfilter:  We suggest  -f "QUAL > 30 " (Depth over 10 reads), "DP > 10" (minimum Phred-scaled probability of error over 30) (Galaxy Version 1.0.0_rc1+galaxy1)
4) VcfAllelicPrimitives: Split allelic primitives into multiple VCF lines, default settings (Galaxy Version 0.0.3)
5) VCF-VCF Intersect: on WT and Null VCF files, mm10 build (Galaxy Version 1.0.0_rc1+galaxy0)


Script Outline

Inputted variants are binned every 10 million base pairs according to its chromosomal coordinates (mm10 build: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26) and ordered in a contingency table. After this, frequency distribution of variants is tested by applying the Cochran-Armitage test for trend distribution, available in the DescTools package in R. (https://cran.r-project.org/web/packages/DescTools/index.html). The program will generate a genome-wide plot of variants per genotype based on the ggplot2 R package (https://cran.r-project.org/web/packages/ggplot2/index.html) and a summary of chromosomes containing KO/KI-ligated variants, based on the frequency distribution of wild-type and KO/KI genotypes. 

R Dependences

The following R/Bioconductor (https://www.bioconductor.org/) packages are required for the script usage

dplyr_0.7.4       
gridExtra_2.3     
reshape2_1.4.2    
ggplot2_2.2.1    
DescTools_0.99.23

Usage
1) Set Directory of VCF files (wild-type and KO/KI, respectively)
2) Open genotype_variants.R file
3) Execute all lines
4) A prompted window will ask for the wild-type VCF file. Select it. 
5) A prompted window will ask for the KO/KI VCF file. Select it
6) Outputs: Genome-wide plot of variants per genotype (PDF file) and a summary of the Cochran-Armitage Test per chromosome (PDF file)

Contributors

genotype_variants.R is made by Carlos Farkas

Support

email

Please report any issues or questions by email to cfarkas@udec.cl or carlosfarkas@gmail.com.

Test 

We provide two Filtered VCF files from a RNA-seq analysis of MEFs (Mouse embryonic fibroblasts) isolated from STC1+/+ and STC1-/- mice https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47395 under GeoDataset accesion GSE47395. 


