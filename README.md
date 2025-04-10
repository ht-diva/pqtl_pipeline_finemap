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
We performed colocalization (Giambartolomei et al., 2014) across the pQTL signals. To meet the fundamental assumption of colocalization of only one causal variant per locus, we used conditional datasets, thus performing one colocalization test per pair of independent SNPs in 2 overlapping loci. For each regional association and each target SNP, we identified a credible set as the set of variants with posterior inclusion probability (PIP) > 0.99 within the region. More precisely, using the conditional dataset, we computed Approximate Bayes Factors (ABF) with the ‘process.dataset’ function in the coloc v5.2.3 R package and calculated posterior probabilities by normalizing ABFs across variants. Variants were ranked, and those with a cumulative posterior probability exceeding 0.99 were included in the credible sets. Among XXX protein pairs with overlapping loci, XXX protein pairs sharing a credible set variant were then tested for colocalization using the ‘coloc.abf’ function. Colocalized pairs were identified when the posterior probability for hypothesis 4 assuming a shared causal variant for two proteins exceeded 0.80.


### New Features on Top of NF pipeline
We also incorporated new features such as exclusion of signals in *HLA* and *NLRP12* regions from the results and follow-up analyses, allowing user to decide through the configuration file.


### NOTE
This SMK pipeline which is designed for pQTLs project **does not** include munging and alignment of input GWAS summary files. Therefore, it is a MUST to have your GWAS results completely harmonized by your genotype data. Eg. variants IDs, refrence/alternate (effect/other) alleles should be concordant across your input files. Our GWAS summary stats from REGENIE are already aligned with QC pipeline (adopted by GWASLab) developed by pQTL analysts team at Health Data Science Center.


### How to run the pipeline:
You can use the default configuration file in config/config_finemap.yaml. Otherwise, prepare your configuration in config/ folder. Then, make sure that `configfile` in workflow/Snakefile matches with your newly created config file name. Then, run the pipeline by typing below command in bash.

```bash
sbatch submit.sh
```

### Not interested to run colocalization?
If you want to **skip running colocalization** with your traits, uncomment this `#--until collect_credible_sets` in Makefile. If you want to skip both COJO and colocalization and only run locus breaker, then change previous option in Makefile to `--until collect_loci` and run the pipeline as mentioned before.

## Workflow example

<img src="dag.svg" alt="example workflow">
