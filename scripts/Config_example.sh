#!/bin/bash

bash ./sort_bam.sh
bash ./variant_collection.sh /path/to/mm10.fa 55
bash ./filtering_combined_mouse.sh 
cd vcf_outputs/
bash ./genotype_variants_mouse.sh WT.intersection.vcf KO.intersection.vcf /path/to/mm10.fa

###
