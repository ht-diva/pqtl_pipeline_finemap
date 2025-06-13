
import pandas as pd
from pathlib import Path

path_coloc = "coloc_res.tsv"
genotype = "/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/bed/qc_recoded_harmonised/impute_recoded_selected_sample_filter_hq_var_new_id_alleles_all"

df = pd.read_csv(path_coloc, sep = "\t")

# Create a new column by concatenating top_cond_a and top_cond_b
df['top_cond_ab'] = df['top_cond_a'] + '_' + df['top_cond_b']


my_df = (
    pd.DataFrame(df, columns=["top_cond_a", "top_cond_b", "top_cond_ab"])
    .drop_duplicates(subset='top_cond_ab')
    .set_index("top_cond_ab", drop=False)
    .sort_index()
)


rule all:
    input:
        expand(str(Path("results/ld_test/{top_snps}.ld")), top_snps=my_df.top_cond_ab),
        "results/ld_test/combined_ld.tsv",
        "results/ld_test/combined_coloc_with_ld.tsv",


rule compute_ld:
    input:
        path_coloc = "coloc_res.tsv"
    output:
        ld  = Path("results/ld_test/{top_snps}.ld"),
    params:
        snp_pair = "{top_snps}",
        bfile = genotype
    conda:
        "envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        source /center/healthds/singularity_functions

        out_ld=$(echo {output.ld} | sed 's/.log//')
        snp_a=$(echo {params.snp_pair} | cut -d'_' -f1)
        snp_b=$(echo {params.snp_pair} | cut -d'_' -f2)

        plink2 \
            --bfile {params.bfile} \
            --ld "$snp_a" "$snp_b" \
            --allow-no-sex \
            --out $out_ld \
            --threads 2  \
            --memory 4000 'require'

        mv {output.ld}.log {output.ld}
        """


# Rule to summarize the results
rule combine_ld:
    input:
        ld  = expand("results/ld_test/{top_snps}.ld", top_snps=my_df.top_cond_ab),
    output:
        summary= "results/ld_test/combined_ld.tsv",
    shell:
        """
        echo -e 'SNP_A\tSNP_B\tr2\tDprime' > {output.summary}
        for ld_file in {input.ld}; do
            pair=$(basename "$ld_file" .ld)
            SNP_A=$(echo $pair | cut -d'_' -f1)
            SNP_B=$(echo $pair | cut -d'_' -f2)
            R2_LINE=$(grep -E 'r\^2\s*=' "$ld_file")
            if [[ -n "$R2_LINE" ]]; then
                R2=$(echo "$R2_LINE" | awk '{{print $3}}')
                Dp=$(echo "$R2_LINE" | awk '{{print $6}}')
            else
                R2="NA"
                Dp="NA"
            fi
            echo -e "$SNP_A\t$SNP_B\t$R2\t$Dp" >> {output.summary}
        done
        """

rule append_ld:
    input:
        coloc = "coloc_res.tsv",
        summary = "results/ld_test/combined_ld.tsv",
    output:
        coloc_ld = "results/ld_test/combined_coloc_with_ld.tsv",
    conda:
        "envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        Rscript workflow/scripts/s10_append_ld.R \
            --coloc {input.coloc} \
            --ld {input.summary} \
            --ofile {output.coloc_ld}
        """
