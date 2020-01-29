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

At this step, users can to execute the genotype_variants.R script (on windows or Mac OSX) to visualize variants in the intersected VCF files.

6) VCFAnnotateGenotypes: on KO intersected VCF files using WT intersected VCF file
7) Filter and Sort: select lines that NOT match "added-genotypes"

Output: KO-linked variants by chromosome 

## R Script Outline:

Inputted variants are binned every 10 million base pairs according to its chromosomal coordinates (mm10 build: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26) and ordered in a contingency table. After this, frequency distribution of variants is tested by applying the Cochran-Armitage test for trend distribution, available in the DescTools package in R. (https://cran.r-project.org/web/packages/DescTools/index.html). The program will generate a genome-wide plot of variants per genotype based on the ggplot2 R package (https://cran.r-project.org/web/packages/ggplot2/index.html) and a summary of chromosomes containing KO/KI-linked variants, based on the frequency distribution of wild-type and KO/KI variants.

## Preeliminars:
### Obtaining and installing R (>=3.5.0)
See https://cloud.r-project.org/ for R installation in linux/ubuntu/windows/(Mac) OS X. R version 3.2.3 comes from default in Ubuntu 16.04 LTS but users with older Ubuntu distributions must upgrade R. A way accomplish this can be the following:
```
# Removing R from system
sudo apt-get remove r-base-core

# Editing sources.list 
sudo su
echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Installing R (version 3.6.0)
sudo apt update; sudo apt install r-base
exit
R
```

The following R packages are required for the script usages:

dplyr >=0.7.4       
gridExtra >=2.3     
reshape2 >=1.4.2    
ggplot2 >=2.2.1
DescTools >=0.99

These packages can be installed in R (macOS/ubuntu) by opening R and typying:
>install.packages("dplyr")<br/>install.packages("gridExtra")<br/>install.packages("reshape2")<br/>install.packages("ggplot2")
<br/>install.packages("DescTools")

## Usage():

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

### Obtaining and installing Freebayes and vcflib:
Clone Freebayes folder in current directory: 
```
git clone --recursive git://github.com/ekg/freebayes.git

# If this line not work, try:

git config --global url.https://github.com/.insteadOf git://github.com/
git clone --recursive git://github.com/ekg/freebayes.git

#Enter Freebayes directory and make:
cd freebayes
make
sudo make install
sudo cp scripts/* /usr/local/bin/

#To check installation, type in terminal:
freebayes
bamleftalign

## Installing vcflib (in freebayes folder):

cd vcflib
make   

#After make, binaries and scripts can be copied in /usr/local/bin with sudo. In vcflib/ directory:

sudo cp scripts/* /usr/local/bin/
sudo cp bin/* /usr/local/bin/

#To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities
```

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 

>sudo apt-get install bedtools

### Obtaining and installing up-to-date SAMtools, bcftools and htslib (version 1.9)
Old samtools version will not work. Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
```
cd htslib-1.9    # and similarly for bcftools and samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools and bcftools...
sudo cp samtools /usr/local/bin/
```
Then in a terminal type
>samtools<br>bcftools

to check 1.9 versions (using htslib v1.9)

### Obtaining and installing BamTools
Complete instructions can be found in README and in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing. Users with privileges can accomplish with sudo: 

>sudo apt install bamtools

### Obtaining and installing tabix (Version: >=1.2.1)

>sudo apt install tabix

### R dependences
As mentioned, several packages needs to be installed in R (≥ 3.3.0). Open a shell, type R and then type:
>install.packages("dplyr")<br/>install.packages("gridExtra")<br/>install.packages("reshape2")<br/>install.packages("ggplot2")<br/>install.packages("DescTools")<br/>

## Execution:
Place the following scripts in a folder:
> sort_bam.sh<br/>variant_collection.sh<br/>filtering_combined_mouse.sh<br/>genotype_variants_mouse.sh<br/>genotype_variants_mouse_linux.R

The reference genome for freebayes:
> mm10.fa (properly indexed with samtools faidx: e.g.: samtools faidx mm10.fa)

and the renamed BAM files from wild-type (WT) and knockout (KO) genotypes with the "WT" and "KO" prefix. e.g.:
> WT1.bam<br/>WT2.bam<br/>KO1.bam<br/>KO2.bam<br/>
#In the case of just one replicate, BAM files must be renamed "WT.bam" or "KO.bam", correspondingly.

Open a terminal and paste the following lines. Using 45 threads for freebayes and assuming replicates per each genotype:

>bash ./sort_bam.sh 45 <br/>bash ./variant_collection.sh /path/to/mm10.fa 45 <br/>bash ./filtering_combined_mouse.sh<br/>cd vcf_outputs/<br/>bash ./genotype_variants_mouse.sh WT.intersection.vcf KO.intersection.vcf /path/to/mm10.fa

These lines can be also executed in BASH (see and edit Config_example.sh in scripts folder) by simply typying:
>bash ./Config_example.sh

### No replicates in one genotype?:
In this example, we will assume replicates only for the KO genotype. In this case, users need to change the "intersection" prefix for "filtered" (e.g. with a single replicate in WT):
```
bash ./sort_bam.sh 45
bash ./variant_collection.sh /path/to/mm10.fa 45
bash ./filtering_combined_mouse.sh 
cd vcf_outputs/
bash ./genotype_variants_mouse.sh WT.filtered.vcf KO.intersection.vcf /path/to/mm10.fa
```

To perform these analysis, we recommend to increase open file limit to 1000000 in the cluster/workstation in use. Please, see "README_ulimit" file for instructions to accomplish this task.

## Test
1) Download WT, KO1, KO2, KO3 BAM alignments (sub-sampĺed datasets from mouse embryonic fibroblasts) from 
https://usegalaxy.org:/u/carlosfarkas/h/test-sall2-ko-rna-seq-gse123168-1 . 

2) Manually rename the four files as WT.bam, KO1.bam, KO2.bam and KO3.bam

3) Place the four renamed BAM files in the Test folder, open a terminal and run the following line:
> bash ./test.sh 40

Step (3) will run the pipeline using 40 threads as specified. 

## Support

Please report any issues or questions by email to cfarkas@udec.cl or carlosfarkas@gmail.com.
