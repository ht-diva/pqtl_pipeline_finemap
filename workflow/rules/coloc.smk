

rule master_file:
    input:
        expand(ws_path("cojo/{seqid}.sentinel"), seqid=analytes.seqid),
    output:
        master = ws_path("coloc/master_coloc.txt"),
    conda:
        "../envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    script:
        "../scripts/s05_collect_info.R"


rule merry_go_round:
    input:
        master = ws_path("coloc/master_coloc.txt"),
    output:
        sentinel = ws_path("coloc/chr{chr_cs}.sentinel"),
    params:
        pairs = ws_path("coloc/chr{chr_cs}_coloc_pairwise_guide_table.tsv"),
        chr = "{chr_cs}",
    conda:
        "../envs/coloc.yml"
    log:
        ws_path("logs/coloc/chr{chr_cs}_coloc_pairwise.log"),
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        Rscript workflow/scripts/s06_find_overlapping_cs.R \
            --coloc_info_table {input.master} \
            --chr_cs {params.chr}  \
            --ofile {params.pairs}
        
        touch {output.sentinel}
        """


rule run_coloc:
    input:
        sentinel = ws_path("coloc/chr{chr_cs}.sentinel"),
    output:
        ofile = ws_path("coloc/chr{chr_cs}_colocalization.table.all.tsv"),
    conda:
        "../envs/fine_mapping.yml"
    params:
        codes = config.get("path_code"),
        pairs = ws_path("coloc/chr{chr_cs}_coloc_pairwise_guide_table.tsv"),
        chr = "{chr_cs}",
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        Rscript workflow/scripts/s07_coloc.R \
            --pipeline_path {params.codes} \
            --coloc_guide_table {params.pairs} \
            --chr_cs {params.chr}  \
            --ofile {output.ofile}
        """
