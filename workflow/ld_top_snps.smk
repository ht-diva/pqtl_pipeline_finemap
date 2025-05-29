
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

#print(my_df)

def get_top_a(wildcards):
    return str(my_df.loc[wildcards, ["top_cond_a", "top_cond_b"]])


rule all:
    input:
        expand(str(Path("results/ld_test/{top_snps}.sentinel")), top_snps=my_df.top_cond_ab)


rule compute_ld:
    input:
        path_coloc = "coloc_res.tsv"
    output:
        sentinel = Path("results/ld_test/{top_snps}.sentinel")
    params:
        #lambda wildcards: get_top_a(wildcards.top_snps)
        ld  = Path("results/ld_test/{top_snps}_ld"),
        snp_pair = "{top_snps}",
        bfile = genotype
    #log:
    #    ws_path("logs/break/{seqid}.log"),
    conda:
        "envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        workflow/pairwise_ld.sh  {params.snp_pair}  {params.bfile}  {params.ld}
        touch {output.sentinel}
        """
