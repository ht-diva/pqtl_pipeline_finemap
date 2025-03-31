
rule run_cojo:
    input:
        gwas=get_sumstats,
        loci=ws_path("break/{seqid}_loci.csv"),
    output:
        sentinel=touch(ws_path("cojo/{seqid}.sentinel")),
        log=ws_path("logs/cojo/{seqid}.log"),
    conda:
        "../envs/fine_mapping.yml"
    params:
        codes=config.get("path_code"),
        geno=config.get("path_geno"),
        g_build=config.get("build"),
        reg_def=config.get("redefine_region"),
        NLRP12=config.get("NLRP12"),
        ofile=ws_path("cojo/{seqid}"),
        ppp=config.get("thresholds").get("ppp"),
        p3=config.get("thresholds").get("p_cojo"),
        p4=config.get("thresholds").get("p_jumper"),
        p5=config.get("thresholds").get("p_signif"),
        p6=config.get("thresholds").get("p_limit"),
        p_label=config.get("labels").get("p_label"),
        chr_label=config.get("labels").get("chr_label"),
        pos_label=config.get("labels").get("pos_label"),
        snpid_label=config.get("labels").get("snpid_label"),
        ea_label=config.get("labels").get("ea_label"),
        oa_label=config.get("labels").get("oa_label"),
        eaf_label=config.get("labels").get("eaf_label"),
        se_label=config.get("labels").get("se_label"),
        beta_label=config.get("labels").get("beta_label"),
        n_label=config.get("labels").get("n_label"),
    log:
        ws_path("logs/cojo/{seqid}.log"),
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
        # Rscript=`ls /conda-envs/*/bin/Rscript`;
        INPUT_FILE={input.loci};

        # read loci and loop over each locus
        while IFS=, read -r col1 col2 col3 col4_onwards; do
            # Assign the values of the first three columns to variables
            chr=$col1
            beg=$col2
            end=$col3

            # Check if chr is empty and exit the loop if it is
            if [[ -z "$beg" ]]; then
                break
            fi

            # Skip the first line then run the script
            if [[ "$chr" != "chr" ]]; then
                Rscript workflow/scripts/s03_cojo_finemapping.R  \
                --pipeline_path   {params.codes}  \
                --dataset_gwas {input.gwas}  \
                --phenotype_id {params.ofile}  \
                --build   {params.g_build}  \
                --lb_bis  {params.reg_def}  \
                --nlrp12  {params.NLRP12}   \
                --chr   "$chr"  \
                --start "$beg"  \
                --end   "$end"  \
                --bfile {params.geno} \
                --cs_thresh {params.ppp}  \
                --p_cojo   {params.p3}  \
                --p_jumper {params.p4}  \
                --p_signif {params.p5}  \
                --p_limit  {params.p6}  \
                --outdir {params.ofile} \
                --plink2_mem {resources.mem_mb}  \
                --plink2_threads {resources.threads} \
                --chr_label {params.chr_label} \
                --pos_label {params.pos_label} \
                --snpid_label {params.snpid_label} \
                --ea_label {params.ea_label} \
                --oa_label {params.oa_label} \
                --eaf_label {params.eaf_label} \
                --se_label {params.se_label} \
                --beta_label {params.beta_label} \
                --n_label {params.n_label} \
                --p_label   {params.p_label} >> {log}
            fi
        done < "$INPUT_FILE"
        #touch {output.sentinel}
        """

rule collect_credible_sets:
    input:
        expand(ws_path("cojo/{seqid}.sentinel"), seqid=analytes.seqid),
        ws_path("break/collected_loci_excluding_mhc.csv")
    output:
        ofile=ws_path("cojo/collected_credible_sets.csv"),
    conda:
        "../envs/locus_breaker.yml"
    resources:
        runtime=lambda wc, attempt: 999 + attempt * 60,
    script:
        "../scripts/s04_collect_credible_sets.R"


