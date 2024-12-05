
rule master_file:
    input:
        ws_path("cojo/{seqid}.sentinel"),
    output:
        sentinel = ws_path("coloc/{seqid}_master_file.sentinel"),
    params:
        info = ws_path("cojo/{seqid}/finemaping/*_coloc_info_table.tsv"),
        ofile = ws_path("coloc/master_coloc.txt"),
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        cat {params.info} >> {params.ofile}
        touch {output.sentinel}
        """


rule merry_go_round:
    input:
        sentinel = expand(ws_path("coloc/{seqid}_master_file.sentinel"), seqid=analytes.seqid),
    output:
        sentinel = ws_path("coloc/chr{chr_cs}.sentinel"),
    params:
        pairs = ws_path("coloc/chr{chr_cs}_coloc_pairwise_guide_table.tsv"),
        info = ws_path("coloc/master_coloc.txt"),
        chr = "{chr_cs}",
    conda:
        "../envs/coloc.yml"
    log:
        ws_path("logs/coloc/chr{chr_cs}_coloc_pairwise.log"),
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        Rscript workflow/scripts/s05_find_overlapping_cs.R \
            --coloc_info_table {params.info} \
            --chr_cs {params.chr}  \
            --ofile {params.pairs}
        
        touch {output.sentinel}
        """
