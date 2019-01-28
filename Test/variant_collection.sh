#!/bin/bash

### Obtaining the Human Reference Genome:
# Download Human Reference Genome mm10: wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
# Download faToTwoBit script, available for linux and macOSX at http://hgdownload.cse.ucsc.edu/admin/exe/.
# Convert 2bit format to fasta: ./twoBitToFa hg19.2bit hg19.fa
# Index fasta file with samtools: samtools faidx hg19.fa

### Obtaining the Mouse Reference Genome:
# Download Mouse Reference Genome mm10: wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit
# Download faToTwoBit script, available for linux and macOSX at http://hgdownload.cse.ucsc.edu/admin/exe/.
# Convert 2bit format to fasta: ./twoBitToFa mm10.2bit mm10.fa
# Index fasta file with samtools: samtools faidx mm10.fa

### Obtaining and installing Freebayes:
## Cloning Freebayes folder in current directory 
# git clone --recursive git://github.com/ekg/freebayes.git
## Enter Freebayes directory and make
# cd freebayes/
# make
## To install to e.g. /usr/local/bin (default), type:
# sudo make install
## To check installation, type in terminal
# freebayes
# bamleftalign

### Obtaining and installing vcflib:
## Cloning vcflib folder in current directory
# git clone --recursive git://github.com/vcflib/vcflib.git
## Enter vcflib directory and make
# cd vcflib/
# make   #If you want to use threading type make openmp instead of make. Only a few VCFLIB tools are threaded.
## After make, binaries and scripts can be copied in /usr/local/bin:
# in vcflib/ directory
# cp scripts/* /usr/local/bin/
# cp bin/* /usr/local/bin/
# To check vcflib scripts, type vcf in terminal followed by TAB and display all posibilities

### Obtaining and Installing BEDTools
## Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html
## Users with privileges can accomplish with sudo: sudo apt-get install bedtools

### Obtaining and installing SAMtools
## Complete instructions can be found in README and in http://www.htslib.org/
## Users with privileges can accomplish with sudo: sudo apt-get install samtools

### Obtaining and installing BamTools
## Complete instructions can be found in README and in https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
## Users with privileges can accomplish with sudo: sudo apt install bamtools 

### Obtaining and installing bcftools
## Complete instructions can be found in README and in https://samtools.github.io/bcftools/
## Users with privileges can accomplish with sudo: sudo apt install bcftools 


REF=${1}
THREADS=${2}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {REFERENCE} {THREADS}"
  echo ""
  echo "This script will call variants using freebayes-parallel in every .sorted.bam file present in the folder"
  echo ""
  echo "REFERENCE: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "The reference genome must be also indexed. e.g.: samtools faidx hg19.fa"
  echo ""
  echo "THREADS: Number of CPUs for the task (integer)"
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {REFERENCE} {THREADS}"
  echo ""
  echo "This script will call variants using freebayes-parallel in every .sorted.bam file present in the folder"
  echo ""
  echo "REFERENCE: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "The reference genome must be also indexed. e.g.: samtools faidx hg19.fa"
  echo ""
  echo "THREADS: Number of CPUs for the task (integer)"
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {REFERENCE} {THREADS}"; exit 1; }

if [ $# -ne 2 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {REFERENCE} {THREADS}"
  exit 3
fi

begin=`date +%s`
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

echo ""
echo "Performing Variant Calling with freebayes for the following bam files:"
echo ""
a= ls -1 *.sorted.bam 
echo ""
echo ${2} "Threads will be used for the task"
echo ""
echo "The output directory will be the following:"
echo ${dir1}
echo ""

for a in *.sorted.bam; do freebayes-parallel <(fasta_generate_regions.py ${REF}.fai 100000) ${THREADS} -f ${REF} -b ${a} > ${a}.vcf
done 
echo ""
end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

### 
