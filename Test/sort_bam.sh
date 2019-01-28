#!/bin/bash

Threads=${1}

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {Threads}"; exit 1; }

for file in *.bam
do
echo $file
bam_name=$(echo ${file} | sed "s/.bam//")
echo $bam_name

# Sort BAM file
samtools sort ${bam_name}.bam > ${bam_name}.sorted.bam --threads ${Threads}

# index the bam file
samtools index ${bam_name}.sorted.bam
done
