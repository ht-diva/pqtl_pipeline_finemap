
rule break_locus:
    input:
        gwas=get_sumstats,
    output:
        loci=ws_path("break/{seqid}_loci.csv"),
    conda:
        "../envs/locus_breaker.yml"
    params:
        codes=config.get("path_code"),
        phenotype_id="{seqid}",
        outdir=ws_path("break/{seqid}"),
        p1=config.get("thresholds").get("p1"),
        p2=config.get("thresholds").get("p2"),
        hole=config.get("thresholds").get("hole"),
        p_label=config.get("labels").get("p_label"),
        chr_label=config.get("labels").get("chr_label"),
        pos_label=config.get("labels").get("pos_label"),
    log:
        ws_path("logs/break/{seqid}.log"),
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    shell:
        """
        # Rscript=`ls /conda-envs/*/bin/Rscript`;
        Rscript workflow/scripts/s01_locus_breaker.R \
            --pipeline_path {params.codes} \
            --input {input.gwas} \
            --phenotype_id '{params.phenotype_id}' \
            --p_thresh1 {params.p1} \
            --p_thresh2 {params.p2} \
            --hole      {params.hole} \
            --outdir    {params.outdir} \
            --p_label   {params.p_label} \
            --chr_label {params.chr_label} \
            --pos_label {params.pos_label} > {log}
        """


rule collect_loci:
    input:
        expand(ws_path("break/{seqid}_loci.csv"), seqid=analytes.seqid),
    output:
        ofile=ws_path("break/collected_loci_excluding_mhc.csv"),
    conda:
        "../envs/locus_breaker.yml"
    params:
        NLRP12=config.get("NLRP12"),
        build=config.get("build")
    resources:
        runtime=lambda wc, attempt: attempt * 30,
    script:
        "../scripts/s02_collect_loci.R"
