# pqtl_pipeline_finemap
Fine mapping analysis within the pQTL pipeline project at Human Technopole, Milan, Italy

We started this analysis pipeline in early April 2024. We adopted the Next-Flow (NF) pipeline developed by the Statistical Genomics team at Human Technopole and deployed it in Snakemake (SMK). We independtly validated each of the multiple analyses stated below before incorporating it in SMK.


### Locus Breaker

We incorporated **Locus Breaker (LB)** function written in R (see publication PMID:) for example meta-analysis GWAS results of the proteins and we deployed it in SMK in mid April 2024.


### COJO Conditional Analysis
Once running the pipeline, `rule run_cojo` will generate output files below:
- list of independent variants resulted from GCTA cojo-slct (TSV/CSV)
- conditional dataset for each independent signal resulted from GCTA cojo-cond (RDS)
- fine-mapping results using coloc::coloc.ABF function, containing values such as l-ABF, posterior probabilities (PPI) for each variant (RDS)
- colocalization info table containing credible set variants (with cumulative PPI > 0.99) for each independent variant
- regional association plots

These outputs are going to be stored in `workspace_path` provided by the user in config_finemap.yaml and stored in such directory: 
<workspace_path>/results/*/cojo/<seqid>


### Colocalization of Two Proteins


### New Features on Top of NF pipeline
We also incorporated new features such as exclusion of signals in *HLA* and *NLRP12* regions from the results and follow-up analyses, allowing user to decide through the configuration file.


### NOTE
The first two steps of the original NF pipeline (mungigg & alignment) are not used in this SMK pipeline for pQTLs. Therefore, it is a MUST to have your GWAS results completely harmonized by your geneotype data. Eg. varinats IDs, refrence/alternate (effect/other) alleles should be concordant across your input files. Our GWAS summary stats from REGENIE are already aligned with QC pipeline (adopted by GWASLab) developed by pQTL analysts team at Health Data Science Center. to do so, you my need a mapping file like here:
`/exchange/healthds/pQTL/results/INTERVAL/qced_sumstats/table.snp_mapping.tsv.gz`

