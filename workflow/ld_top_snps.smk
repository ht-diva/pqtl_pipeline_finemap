
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
