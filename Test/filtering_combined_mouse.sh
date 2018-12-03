#!/bin/bash

dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "The output directory will be the following:"
echo ""
echo ${dir1}
echo ""
### 
echo "Filtering with QUAL > 30 the following VCF files:"
a= ls -1 *.vcf
for a in *.vcf; do vcffilter -f "QUAL > 30"  ${a} > ${a}_QUAL.vcf
done
###
mkdir vcf_inputs
mv *.sorted.bam.vcf ${dir1}/vcf_inputs/
echo ""
echo ""
echo "Done. Filtering with DP > 10:"
b= ls -1 *_QUAL.vcf
for b in *_QUAL.vcf; do vcffilter -f "DP > 10"  ${b} > ${b}_DP.vcf
done
rm *.sorted.bam.vcf_QUAL.vcf
###
for file in *.sorted.bam.vcf_QUAL.vcf_DP.vcf
do
  mv "$file" "${file/.sorted.bam.vcf_QUAL.vcf_DP.vcf/_filtered.vcf}"
done
###
echo ""
echo "Done. Decomposing the following files:"
c= ls -1 *_filtered.vcf
for c in *_filtered.vcf; do vcfallelicprimitives -g ${c} > ${c}_decomposed.vcf
done
rm *_filtered.vcf
###
for file in *_filtered.vcf_decomposed.vcf
do
  mv "$file" "${file/_filtered.vcf_decomposed.vcf/.filtered.vcf}"
done
echo ""
echo ""
mkdir vcf_outputs
mv *.filtered.vcf ${dir1}/vcf_outputs/
echo "Done. The filtered and decomposed VCF files are in the vcf_outputs folder"
cp genotype_variants_mouse.sh genotype_variants_mouse_linux.R ${dir1}/vcf_outputs/
###
cd ${dir1}/vcf_outputs
echo ""

if [ -f WT.filtered.vcf ]; then
    echo "WT.filtered.vcf file found. Continue with KO vcf files"
    echo ""
    : 
else
    echo "intersecting the following WT VCF files:"
    d= ls -1 WT*.filtered.vcf
    for d in WT*.filtered.vcf; do bgzip -i ${d}
    done
    ###
    echo ""
    echo "BGZIP compressed files:"
    e= ls -1 WT*.filtered.vcf.gz
    for e in WT*.filtered.vcf.gz; do zcat ${e} | bgzip -c > ${e}.new.gz && tabix -f ${e}.new.gz
    done
    echo ""
    echo "Intersecting with bcftools isec:"
    bcftools isec -n =$(ls -1 WT*.filtered.vcf.gz.new.gz | wc -l ) WT*.filtered.vcf.gz.new.gz > WT.sites.vcf
    bgzip -d WT1.filtered.vcf.gz
    ###
    grep "#" WT1.filtered.vcf > header.vcf
    awk 'FILENAME == "WT.sites.vcf" { remember[$1 $2]=1 ;}
    FILENAME != "WT.sites.vcf" { if ( $1 $2 in remember ) print ; } ' WT.sites.vcf WT1.filtered.vcf > WT_1.vcf
    cat header.vcf WT_1.vcf > WT.intersection.vcf
    bgzip -d WT*.filtered.vcf.gz
    rm *.new.gz *.gzi header.vcf WT_1.vcf WT.sites.vcf *.tbi
    echo "Intersection done. The Wild-type intersected file is called WT.intersection.vcf"  
fi
###
###
###
if [ -f KO.filtered.vcf ]; then
    echo "KO.filtered.vcf file found. All done."
    : 
else
    echo "Intersecting the following KO VCF files:"
    d= ls -1 KO*.filtered.vcf
    for d in KO*.filtered.vcf; do bgzip -i ${d}
    done
    ###
    echo ""
    echo "BGZIP compressed files:"
    e= ls -1 KO*.filtered.vcf.gz
    for e in KO*.filtered.vcf.gz; do zcat ${e} | bgzip -c > ${e}.new.gz && tabix -f ${e}.new.gz
    done
    echo ""
    echo "Intersecting with bcftools isec:"
    bcftools isec -n =$(ls -1 KO*.filtered.vcf.gz.new.gz | wc -l ) KO*.filtered.vcf.gz.new.gz > KO.sites.vcf
    bgzip -d KO1.filtered.vcf.gz
    ###
    grep "#" KO1.filtered.vcf > header.vcf
    awk 'FILENAME == "KO.sites.vcf" { remember[$1 $2]=1 ;}
    FILENAME != "KO.sites.vcf" { if ( $1 $2 in remember ) print ; } ' KO.sites.vcf KO1.filtered.vcf > KO_1.vcf
    cat header.vcf KO_1.vcf > KO.intersection.vcf
    bgzip -d KO*.filtered.vcf.gz
    rm *.new.gz *.gzi header.vcf KO_1.vcf KO.sites.vcf *.tbi
    echo "Intersection done. The intersected file is called KO.intersection.vcf"  
    echo "All done. The intersected files are located in the vcf_outputs/ folder"
fi
