#!/usr/bin/bash

snp_pair=$1
bfile=$2
ld=$3

source /center/healthds/singularity_functions

out_ld=$(echo $ld | sed 's/\:/\_/g')
snp_a=$(echo $snp_pair | cut -d'_' -f1)
snp_b=$(echo $snp_pair | cut -d'_' -f2)

plink2 \
    --bfile $bfile \
    --ld "$snp_a" "$snp_b" \
    --allow-no-sex \
    --out $out_ld \
    --threads 2  \
    --memory 4000 'require'
