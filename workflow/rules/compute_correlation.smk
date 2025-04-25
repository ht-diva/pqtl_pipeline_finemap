
from pathlib import Path
import pandas as pd

path_merged = "conf/path_merged.txt"
df = pd.read_csv(path_merged, sep = ",")

my_df = (
    pd.DataFrame(df, columns=["seqid", "path_3m", "path_bb"])
    .set_index("seqid", drop=False)
    .sort_index()
)


def get_path_3m(wildcards):
    return my_df.loc[wildcards, "path_3m"]

def get_path_bb(wildcards):
    return my_df.loc[wildcards, "path_bb"]


rule all:
    input:
        expand(str(Path("results/cor_beta/single/{seqid}.tsv")), seqid = df.seqid),
        "results/cor_beta/combined.tsv"


rule compute_r:
    input:
        script = "workflow/scripts/compute_r.R",
        path_3m = lambda wildcards: get_path_3m(wildcards.seqid),
        path_bb = lambda wildcards: get_path_bb(wildcards.seqid),
    output:
        ofile = Path("results/cor_beta/single/{seqid}.tsv")
    conda:
        "../envs/locus_breaker.yml"
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
        mem_mb=lambda wc, attempt: 3000 + attempt * 1000,
    shell:
        """
        Rscript {input.script} {input.path_3m} {input.path_bb} {output.ofile}
        """


rule combine_rs:
    input:
        expand("results/cor_beta/single/{seqid}.tsv", seqid=df.seqid)
    output:
        "results/cor_beta/combined.tsv"
    run:
        with open(output[0], 'w') as out_file:
            for idx, f in enumerate(input):
                with open(f) as in_file:
                    lines = in_file.readlines()
                    if idx == 0:
                        out_file.writelines(lines)
                    else:
                        out_file.writelines(lines[1:])  # Skip header
