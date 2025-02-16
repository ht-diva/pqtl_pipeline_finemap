#!/bin/bash

vcf_raw=$1
txt_out=$2
chain_file=$3
input_hg37=$4
input_hg38=$5

vcf_zip=$(echo $vcf_raw | sed -E 's/vcf/vcf.gz/')
vcf_std=$(echo $vcf_raw | sed -E 's/vcf/std.vcf.gz/')
vcf_out=$(echo $vcf_raw | sed -E 's/txt/vcf/')

source /exchange/healthds/singularity_functions


# LiftOver the VCF file
bgzip -c ${vcf_raw} > ${vcf_zip} && \
tabix -p vcf ${vcf_zip} && \
bcftools norm -f ${input_hg37} -c s -Oz -o ${vcf_std} ${vcf_zip} && \
bcftools +liftover --no-version -Ou ${vcf_std} -- -s ${input_hg37} -f ${input_hg38} -c ${chain_file} > ${vcf_out} && \
bcftools query -f '%CHROM\t%POS\t%ID' ${vcf_out} > ${txt_out} && \
rm ${vcf_zip} ${vcf_std}
echo "VCF file conversion and liftOver completed successfully!"
