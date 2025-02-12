#!/bin/bash

chain_file=/group/diangelantonio/public_data/liftOver/hg19ToHg38.over.chain.gz
input_hg37=/group/diangelantonio/public_data/liftOver/human_g1k_v37.fasta
input_hg38=/group/diangelantonio/public_data/liftOver/hg38.fa
vcf_raw=$1
vcf_zip=$(echo $vcf_raw | sed -E 's/vcf/vcf.gz/')
vcf_std=$(echo $vcf_raw | sed -E 's/vcf/std.vcf.gz/')
vcf_b38=$(echo $vcf_raw | sed -E 's/vcf/hg38.vcf.gz/')
out_txt=$(echo $vcf_raw | sed -E 's/vcf/txt/')

source /exchange/healthds/singularity_functions

# echo $vcf_raw
# echo $vcf_zip
# echo $vcf_std
# echo $vcf_b38
# echo $out_txt

# LiftOver the VCF file
bgzip -c ${vcf_raw} > ${vcf_zip} && \
tabix -p vcf ${vcf_zip} && \
bcftools norm -f ${input_hg37} -c s -Oz -o ${vcf_std} ${vcf_zip} && \
bcftools +liftover --no-version -Ou ${vcf_std} -- -s ${input_hg37} -f ${input_hg38} -c ${chain_file} > ${vcf_b38} && \
bcftools query -f '%CHROM\t%POS\t%ID' ${vcf_b38} > ${out_txt}
# zcat seq.12571.14_10_104214048_105206365_target_10:104436641:C:T.vcf.gz | grep -v '^##'
