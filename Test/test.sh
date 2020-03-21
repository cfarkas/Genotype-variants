#!/bin/bash

CORES=${1}

# Obtaining mm10.fa genome, and indexing

wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit
chmod 755 twoBitToFa
./twoBitToFa mm10.2bit mm10.fa
samtools faidx mm10.fa


# Pipeline
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
bash ./sort_bam.sh ${CORES}
bash ./variant_collection.sh mm10.fa ${CORES} 
bash ./filtering_combined_mouse.sh 
cd vcf_outputs/
bash ./genotype_variants_mouse.sh WT.filtered.vcf KO.intersection.vcf /${dir1}/mm10.fa

###
