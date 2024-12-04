
rule make_file:
    input:
        expand(ws_path("cojo/{seqid}.sentinel"), seqid=analytes.seqid),
    output:
        ofile = ws_path("coloc/master_coloc.txt"),
    params:
        outdir = config.get("workspace_path")
    log:
        ws_path("logs/coloc/master_file.log"),
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    shell:
        """
        cat {params.outdir}/cojo/{seqid}/finemaping/*_coloc_info_table.tsv >> {output.ofile}
        """

