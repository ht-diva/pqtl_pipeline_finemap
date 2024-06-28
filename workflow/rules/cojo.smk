
rule run_cojo:
    input:
        gwas=get_sumstats,
        loci=ws_path("break/{seqid}_loci.csv"),
    output:
        sentinel=ws_path("cojo/{seqid}/sentinel.txt"),
        log=ws_path("logs/cojo/{seqid}.log"),
    conda:
        "../envs/r_environment.yml"
    params:
        codes=config.get("path_code"),
        geno=config.get("path_geno"),
        ofile=ws_path("cojo/{seqid}"),
        ppp=config.get("thresholds").get("ppp"),
        p3=config.get("thresholds").get("p3"),
        p4=config.get("thresholds").get("p4"),
        maf=config.get("thresholds").get("maf"),
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
        runtime=lambda wc, attempt: attempt * 60,
        mem_mb=lambda wc, attempt: 8000 + attempt * 2048
    shell:
        """
        Rscript=`ls /conda-envs/*/bin/Rscript`;
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
                Rscript ../scripts/s03_cojo_finemapping.R  \
                --pipeline_path   {params.codes}  \
                --dataset_gwas {input.gwas}  \
                --phenotype_id {params.ofile}  \
                --chr   "$chr"  \
                --start "$beg"  \
                --end   "$end"  \
                --bfile {params.geno} \
                --cs_thresh {params.ppp}  \
                --p_thresh3 {params.p3}  \
                --p_thresh4 {params.p4}  \
                --maf       {params.maf}  \
                --outdir {params.ofile} \
                --plink2_mem {resources.mem_mb}  \
                --plink2_threads {threads} \
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
        touch {output.sentinel}
        """
