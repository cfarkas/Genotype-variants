#!/bin/bash

CORES=${1}

# Obtaining mm10.fa genome, and indexing

wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit
./twoBitToFa mm10.2bit mm10.fa
samtools faidx mm10.fa


# Pipeline

bash ./sort_bam.sh
bash ./variant_collection.sh /datos1/genotype_variants_mouse/mm10.fa ${CORES} 
bash ./filtering_combined_mouse.sh 
cd vcf_outputs/
bash ./genotype_variants_mouse.sh WT.filtered.vcf KO.intersection.vcf /datos1/genotype_variants_mouse/mm10.fa

###
