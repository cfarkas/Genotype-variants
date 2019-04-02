# Genotype-variants: A streamlined pipeline for the genetic characterization of Genetically Engineered Mice (GEM) models based on NGS data. 

## Pipeline Outline

Genotype Variants is an automated BASH/R pipeline for the adequate visualization of RNA-seq variants on the mouse mm10 build, based on VCF files (Variant Calling Format files, http://samtools.github.io/hts-specs/VCFv4.3.pdf). The pipeline takes as input BAM files from RNA-Seq, whole exome sequencing (WES) or whole genome sequencing (WGS) from wild-type and KO/KI genotypes. The pipeline call variants in each sampĺe with the Freebayes variant caller, intersect variants from each genotype and outputs filtered VCF files per genotype, genome-wide plots of variants and a summary of chromosomes with differential distribution of variants among genotypes. This pipeline can be fully implemented using the main galaxy tools (https://usegalaxy.org/). 

The input could be BAM files from a splice-aware aligner: e.g. HISAT2 or BAM files from Bowtie2 or BWA (form WES/WGS data)
                  
To collect variants, BAM files from the alignment of RNA-seq/WES/WGS reads are sorted and genome-wide simple diploid calling is performed with Freebayes (see https://github.com/ekg/freebayes). Raw VCF files are further processed using the VCFlib toolkit (https://github.com/vcflib/vcflib). VCF files needs to be intersected among biological replicates from each genotype (VCF-VCF intersect tool) and the resulting VCF file is then filtered using an adequate criteria. Finally, a variant normalization step is performed with the VcfAllelicPrimitives tool from VCFlib, simplifying multi-nucleotides variants (MNPs) into primitive alleles.  
The pipeline can be executed in two ways:

A) Galaxy to obtain the intersected VCF files per genotype and then in R (from windows and Max OSX)
B) Fully automated in BASH by using a simple config file.


# Galaxy Pipeline: 

Go to https://usegalaxy.org/ and create an account. Users can also install Galaxy locally following these instructions: https://galaxyproject.org/admin/get-galaxy/

1) Input: Sorted BAM file from alignment (e.g HISAT)
2) FreeBayes: Simple diploid calling, default settings 
3) VCFfilter:  We suggest  -f "QUAL > 30 " (Depth over 10 reads), "DP > 10" (minimum Phred-scaled probability of error over 30) 
4) VcfAllelicPrimitives: Split allelic primitives into multiple VCF lines, default settings (Galaxy Version 0.0.3)
5) VCF-VCF Intersect: on WT and KO VCF files, mm10 build

At this step, users needs to execute the genotype_variants.R script (on windows or Mac OSX)

## Script Outline:

Inputted variants are binned every 10 million base pairs according to its chromosomal coordinates (mm10 build: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26) and ordered in a contingency table. After this, frequency distribution of variants is tested by applying the Cochran-Armitage test for trend distribution, available in the DescTools package in R. (https://cran.r-project.org/web/packages/DescTools/index.html). The program will generate a genome-wide plot of variants per genotype based on the ggplot2 R package (https://cran.r-project.org/web/packages/ggplot2/index.html) and a summary of chromosomes containing KO/KI-linked variants, based on the frequency distribution of wild-type and KO/KI variants.

## Preeliminars (R dependences):

The following R packages are required for the script usage (see https://www.r-project.org/ for R installation)

dplyr >=0.7.4       
gridExtra >=2.3     
reshape2 >=1.4.2    
ggplot2 >=2.2.1    
DescTools >=0.99.23

These packages can be installed in R (windows/macOS/ubuntu) by opening R and typying:
>install.packages("dplyr")<br/>install.packages("gridExtra")<br/>install.packages("reshape2")<br/>install.packages("ggplot2")<br/>install.packages("DescTools")<br/>

To check installation of these packages, open R and type:
>library(dplyr)<br/>library(gridExtra)<br/>library(reshape2)<br/>library(ggplot2)<br/>library(DescTools)<br/>

## Usage:

1) Place genotype_variants.R and the VCF files (wild-type and KO/KI outputs after the galaxy pipeline) in the same folder<br/>
2) Open genotype_variants.R file through R (go to File --> Open Script/Document)<br/>
3) Execute all lines (Edit --> Select all, and then Edit --> Execute)<br/>
4) A prompted window will ask for the wild-type VCF file. Select it.<br/>
5) A prompted window will ask for the KO/KI VCF file. Select it.<br/>
6) Output: Genome-wide plot of variants per genotype (PDF file) and a summary of the Cochran-Armitage Test per chromosome (PDF file)



# BASH pipeline:

Users can execute the pipeline by using several scripts sequentially:<br/>

sort_bam.sh (for sorting BAM files)<br/>variant_collection.sh (for calling variants with freebayes)<br/>filtering_combined_mouse.sh (filtering and intersection of VCF files per genotype)<br/>genotype_variants_mouse.sh (KO-ligated variants and R plots generation).<br/>

These four scripts should be present in the working folder along with the genotype_variants_mouse_linux.R script, the mm10.fa genome (in FASTA format, properly indexed) and the BAM files (unsorted) from every genotype to analyze. 

## Preeliminars:

### Obtaining the Mouse Reference Genome:
Download Mouse Reference Genome mm10: 
>wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit

Download faToTwoBit script, available for linux and macOSX at http://hgdownload.cse.ucsc.edu/admin/exe/<br/>
Users needs to change script permissions and convert 2bit format to fasta:
> chmod 755 twoBitToFa<br/>
>./twoBitToFa mm10.2bit mm10.fa

Finally, index fasta file with samtools: 
>samtools faidx mm10.fa

### Obtaining and installing Freebayes:
Clone Freebayes folder in current directory: 
>git clone --recursive git://github.com/ekg/freebayes.git

Enter Freebayes directory and make:
>cd freebayes/<br/>make

To install to e.g. /usr/local/bin (default), type:
>sudo make install

To check installation, type in terminal:
>freebayes<br/>bamleftalign

### Obtaining and installing vcflib:
Clone vcflib folder in current directory:
>git clone --recursive git://github.com/vcflib/vcflib.git

Enter vcflib directory and make
>cd vcflib/<br/>make   

After make, binaries and scripts can be copied in /usr/local/bin with sudo. In vcflib/ directory:

>sudo cp scripts/* /usr/local/bin/<br/>sudo cp bin/* /usr/local/bin/

To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 

>sudo apt-get install bedtools

### Obtaining and installing SAMtools
Complete instructions can be found in README and in http://www.htslib.org/. Users with privileges can accomplish with sudo: 

>sudo apt-get install samtools

### Obtaining and installing BamTools
Complete instructions can be found in README and in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing. Users with privileges can accomplish with sudo: 

>sudo apt install bamtools 

### Obtaining and installing bcftools
Complete instructions can be found in README and in https://samtools.github.io/bcftools/. Users with privileges can accomplish with sudo: 

>sudo apt install bcftools

### R dependences
As mentioned, several packages needs to be installed in R. Open a shell, type R and then type:
>install.packages("dplyr")<br/>install.packages("gridExtra")<br/>install.packages("reshape2")<br/>install.packages("ggplot2")<br/>install.packages("DescTools")<br/>

## Execution:
Place the following scripts in a folder:
> sort_bam.sh<br/>variant_collection.sh<br/>filtering_combined_mouse.sh<br/>genotype_variants_mouse.sh<br/>genotype_variants_mouse_linux.R

The reference genome for freebayes:
> mm10.fa (indexed with samtools faidx)

and the renamed BAM files from wild-type (WT) and knockout (KO) genotypes with the "WT" and "KO" prefix. e.g.:
> WT1.bam<br/>WT2.bam<br/>KO1.bam<br/>KO2.bam<br/>
#In the case of just one replicate, BAM files must be renamed "WT.bam" or "KO.bam", correspondingly.

Open a terminal and paste the following lines. Using 45 threads for freebayes and assuming replicates per each genotype:

>bash ./sort_bam.sh 45 <br/>bash ./variant_collection.sh /path/to/mm10.fa 45 <br/>bash ./filtering_combined_mouse.sh<br/>cd vcf_outputs/<br/>bash ./genotype_variants_mouse.sh WT.intersection.vcf KO.intersection.vcf /path/to/mm10.fa

Using 45 threads for freebayes and assuming replicates only for the KO genotype. In this case, users need to change the "intersection" prefix for "filtered" (e.g. with a single replicate in WT):

>bash ./sort_bam.sh 45 <br/>bash ./variant_collection.sh /path/to/mm10.fa 45<br/>bash ./filtering_combined_mouse.sh <br/>cd vcf_outputs/<br/>bash ./genotype_variants_mouse.sh WT.filtered.vcf KO.intersection.vcf /path/to/mm10.fa

These lines can be also executed in BASH (see and edit Config_example.sh in scripts folder) by simply typying:
>bash ./Config_example.sh

To perform these analysis, we recommend to increase open file limit to 1000000 in the cluster/workstation in use. Please, see "README_ulimit" file for instructions to accomplish this task.

## Test
1) Download WT, KO1, KO2, KO3 BAM alignments (sub-sampĺed datasets from mouse embryonic fibroblasts) from 
https://usegalaxy.org:/u/carlosfarkas/h/test-sall2-ko-rna-seq-gse123168-1 . 

2) Manually rename the four files as WT.bam, KO1.bam, KO2.bam and KO3.bam

3) Place the four renamed BAM files in the Test folder, open a terminal and run the following line:
> bash ./Config_test.sh 40

Step (3) will run the pipeline using 40 threads as specified. 

## Support

Please report any issues or questions by email to cfarkas@udec.cl or carlosfarkas@gmail.com.
