#!/bin/bash

WT_variants=${1}
KO_variants=${2}
REF=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {WT intersected variants} {KO intersected variants} {REFERENCE}"
  echo ""
  echo "This script will invoke an R script to plot variants from WT and KO intersected vcf files"
  echo ""
  echo "WT intersected variants: WT.intersection.vcf file (multiple replicates) or WT.filtered.vcf file (single replicate)"
  echo ""
  echo "KO intersected variants: KO.intersection.vcf file (multiple replicates) of KO.filtered.vcf file (single replicate)"
  echo ""
  echo "REFERENCE: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: bash ./`basename $0` {WT intersected variants} {KO intersected variants} {REFERENCE}"
  echo ""
  echo "This script will invoke and R script to plot variants from WT an KO intersected vcf files"
  echo ""
  echo "WT intersected variants: WT.intersection.vcf file (multiple replicates) or WT.filtered.vcf file (single replicate)"
  echo ""
  echo "KO intersected variants: KO.intersection.vcf file (multiple replicates) of KO.filtered.vcf file (single replicate)"
  echo ""
  echo "REFERENCE: PATH where the reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: bash ./`basename $0` {WT intersected variants} {KO intersected variants} {REFERENCE}"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: bash ./`basename $0` {WT intersected variants} {KO intersected variants} {REFERENCE}"
  exit 3
fi

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
mkdir R_inputs
cp ${1} ${dir1}/R_inputs/
cp ${2} ${dir1}/R_inputs/
mv "${1}" "${1/*.vcf/WT_variants}"
mv "${2}" "${2/*.vcf/KO_variants}"
Rscript genotype_variants_mouse_linux.R 
rm WT_variants
rm KO_variants
mkdir R_outputs
cp graph.pdf ${dir1}/R_outputs/
cp Summary.pdf ${dir1}/R_outputs/
rm graph.pdf Summary.pdf Rplots.pdf
cd ${dir1}/R_inputs/
vcfintersect -i ${1} ${2} -r ${REF} --invert > KO-ligated-variants.vcf
mkdir ${dir1}/KO_ligated_variants
cp ${dir1}/R_inputs/KO-ligated-variants.vcf ${dir1}/KO_ligated_variants/KO-ligated-variants.vcf
rm ${dir1}/R_inputs/KO-ligated-variants.vcf
echo ""
echo "The KO-ligated variants are located in the /vcf_outputs/KO_ligated_variants/ directory"



