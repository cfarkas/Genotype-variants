# Genotype-variants: An streamlined pipeline for the genetic characterization of Genetically Engineered Mice (GEM) models based on NGS data. 

Pipeline Outline

Genotype Variants is an automated BASH/R pipeline for the adequate visualization of RNA-seq variants on the mouse mm10 build, based on VCF files (Variant Calling Format files, http://samtools.github.io/hts-specs/VCFv4.3.pdf). The pipeline takes as input BAM files from RNA-Seq, whole exome sequencing (WES) or whole genome sequencing (WGS) from wild-type and KO/KI genotypes. The pipeline call variants in each sampĺe with the Freebayes variant caller, intersect variants from each genotype and outputs filtered VCF files per genotype, genome-wide plots of variants and a summary of chromosomes with differential distribution of variants among genotypes. This pipeline can be fully implemented using the main galaxy tools (https://usegalaxy.org/). 

The input could be BAM files from a splice-aware aligner: e.g. HISAT2 or BAM files from Bowtie2 or BWA (form WES/WGS data)
                  
To collect variants, BAM files from the alignment of RNA-seq/WES/WGS reads are sorted and genome-wide simple diploid calling is performed with Freebayes (see https://github.com/ekg/freebayes). Raw VCF files are further processed using the VCFlib toolkit (https://github.com/vcflib/vcflib). VCF files needs to be intersected among biological replicates from each genotype (VCF-VCF intersect tool) and the resulting VCF file is then filtered using an adequate criteria. Finally, a variant normalization step is performed with the VcfAllelicPrimitives tool from VCFlib, simplifying multi-nucleotides variants (MNPs) into primitive alleles.  
The pipeline can be executed in two ways:

A) Galaxy to obtain the intersected VCF files per genotype and then in R (from windows and Max OSX)
B) Fully automated in BASH by using a simple config file.

##################################################################################################
##################################################################################################

# Galaxy Pipeline: 

Go to https://usegalaxy.org/ and create an account. Users can also install Galaxy locally following these instructions: https://galaxyproject.org/admin/get-galaxy/

1) Input: Sorted BAM file from alignment (e.g HISAT)
2) FreeBayes: Simple diploid calling, default settings 
3) VCFfilter:  We suggest  -f "QUAL > 30 " (Depth over 10 reads), "DP > 10" (minimum Phred-scaled probability of error over 30) 
4) VcfAllelicPrimitives: Split allelic primitives into multiple VCF lines, default settings (Galaxy Version 0.0.3)
5) VCF-VCF Intersect: on WT and Null VCF files, mm10 build

At this step, users needs to execute the genotype_variants.R script (on windows or Mac OSX)

#Script Outline:

Inputted variants are binned every 10 million base pairs according to its chromosomal coordinates (mm10 build: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26) and ordered in a contingency table. After this, frequency distribution of variants is tested by applying the Cochran-Armitage test for trend distribution, available in the DescTools package in R. (https://cran.r-project.org/web/packages/DescTools/index.html). The program will generate a genome-wide plot of variants per genotype based on the ggplot2 R package (https://cran.r-project.org/web/packages/ggplot2/index.html) and a summary of chromosomes containing KO/KI-ligated variants, based on the frequency distribution of wild-type and KO/KI genotypes. 

# R dependences:

The following R/Bioconductor packages are required for the script usage (see https://www.bioconductor.org/ for installation)

dplyr >=0.7.4       
gridExtra >=2.3     
reshape2 >=1.4.2    
ggplot2 >=2.2.1    
DescTools >=0.99.23

# Usage:

1) Place genotype_variants.R and the VCF files (wild-type and KO/KI outputs after the galaxy pipeline) in the same folder
2) Open genotype_variants.R file through R (go to File --> Open Script/Document)
3) Execute all lines (Edit --> Select all, and then Edit --> Execute)
4) A prompted window will ask for the wild-type VCF file. Select it. 
5) A prompted window will ask for the KO/KI VCF file. Select it
6) Output: Genome-wide plot of variants per genotype (PDF file) and a summary of the Cochran-Armitage Test per chromosome (PDF file)

##################################################################################################
##################################################################################################

### BASH pipeline:

Users can execute the pipeline in bash following a config file. This file will invoke four scripts (sort_bam.sh, variant_collection.sh, filtering_combined_mouse.sh and genotype_variants_mouse.sh) that should be present in the working folder along with the genotype_variants_mouse_linux.R script, the mm10.fa genome (properly indexed) and the BAM files (unsorted) from every genotype to analyze.  

### Preeliminars:

## Obtaining the Mouse Reference Genome:
# Download Mouse Reference Genome mm10: wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit
# Download faToTwoBit script, available for linux and macOSX at http://hgdownload.cse.ucsc.edu/admin/exe/.
# Convert 2bit format to fasta: ./twoBitToFa mm10.2bit mm10.fa
# Index fasta file with samtools: samtools faidx mm10.fa

## Obtaining and installing Freebayes:
## Cloning Freebayes folder in current directory 
# git clone --recursive git://github.com/ekg/freebayes.git
## Enter Freebayes directory and make
# cd freebayes/
# make
## To install to e.g. /usr/local/bin (default), type:
# sudo make install
## To check installation, type in terminal:
# freebayes
# bamleftalign

## Obtaining and installing vcflib:
## Cloning vcflib folder in current directory
# git clone --recursive git://github.com/vcflib/vcflib.git
## Enter vcflib directory and make
# cd vcflib/
# make   #If you want to use threading type make openmp instead of make. Only a few VCFLIB tools are threaded.
## After make, binaries and scripts can be copied in /usr/local/bin with sudo:
# in vcflib/ directory
# sudo cp scripts/* /usr/local/bin/
# sudo cp bin/* /usr/local/bin/
# To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities

## Obtaining and Installing BEDTools
# Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html
# Users with privileges can accomplish with sudo: sudo apt-get install bedtools

## Obtaining and installing SAMtools
# Complete instructions can be found in README and in http://www.htslib.org/
# Users with privileges can accomplish with sudo: sudo apt-get install samtools

## Obtaining and installing BamTools
# Complete instructions can be found in README and in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
# Users with privileges can accomplish with sudo: sudo apt install bamtools 

## Obtaining and installing bcftools
# Complete instructions can be found in README and in https://samtools.github.io/bcftools/
# Users with privileges can accomplish with sudo: sudo apt install bcftools 

### Execution:


Contributors

genotype_variants.R is made by Carlos Farkas

Support


Please report any issues or questions by email to cfarkas@udec.cl or carlosfarkas@gmail.com.

Test 

We provide two Filtered VCF files from a RNA-seq analysis of MEFs (Mouse embryonic fibroblasts) isolated from STC1+/+ and STC1-/- mice https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47395 under GeoDataset accesion GSE47395. 


