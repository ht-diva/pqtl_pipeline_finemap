# pqtl_pipeline_finemap
Fine mapping analysis within the pQTL pipeline project

```{bash}

# get to the project dir
cd projects/pqtl_pipeline_finemap/

# copy R scripts into home directory
cp /exchange/healthds/pQTL/pQTL_workplace/fine_mapping_genomics/nf-hcoloc/bin/*.R projects/pqtl_pipeline_finemap/workflow/scripts/

# create conda env to install and load R libraries
mkdir workflow/envs
cd envs/
touch r_finemap.yml

# add packages with their version in yaml file
conda env create -f r_environment.yml

# test r code
conda activate r_finemap
Rscript ../scripts/s01_sumstat_munging.R
conda deactivate

# get to the env to install missed libraries
cd ../envs/
conda env update -f r_environment.yml
conda env list

# try to reactivate env to test it
conda activate r_finemap
Rscript ../scripts/s01_sumstat_munging.R
conda deactivate

```

- First step is tested. Second is ongoing! (19:00, Fri, 12-Apr-24)

- Configuration and snakefile were created to integrate the codes into the pipeline (23:50, Sun, 14-Apr-24)!

- Lunching the three first steps of the workflow on slurm integrated into snakemake pipeline (03:15, Mon, 15-Apr-24).

- Step four was also debugged after integrating to the pipeline (16:30, Mon, 22-Apr-24).

- It was critical and needed to loop the fine-mapping analysis (step 4) over the identified index variants (P<5e-8, step 3 or locus breaker output) for each protein sequence. This action was done inside rule 'run_cojo' which allowed iterating the R script `s04_cojo_finemapping.R` inside the rule rather than looping the rule through the index variants. Therefore, we read TSV file containing the index variants and iterating the R script only for those variants. So that if there is no genome-wide significant variants for a protein sequence, rule 'run_cojo' would fail for that proetin.

- According to the abovementioned explanations, fine-mapping outputs including independet variants, credible set, regional association plots were saved for 3 out of 6 tested protein sequence in `~/projects/pqtl_pipeline_finemap/output/cojo/<seqid>` (Thu, 23:45, 25-Apr-24).

- The first two steps of the original workflow are not needed for the pQTLs pipeline. Therefore, it is better to start with the step 3, locus breaker, directly using protein GWAS summary stats from REGENIE which is already aligned (Fri, 11:58, 26-Apr-24).

- Locus breaker is incorporated into the pipeline. Time to do so with COJO execution (Sat, 03:45, 27-Apr-24).

- Testing the pipeline on the meta-analysis results modified by GWASLab (Sat, 16:46, 11-May-24).
```bash
ls -1 /exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats/output/seq*/seq*.gwaslab.tsv.gz > conf/path_list_meta.txt
```

```bash

# mapping files
head mapping_chr6.txt
zcat /exchange/healthds/pQTL/results/INTERVAL/qced_sumstats/table.snp_mapping.tsv.gz

# Gcount
cat ThicxTShY15uHN8Y6S44_sum.txt | grep 6:32561656

# Metal
zcat /exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats/output/seq.13435.31/seq.13435.31.gwaslab.tsv.gz

# subset genotype
cat 6_31305620_32886291/ThicxTShY15uHN8Y6S44.bim | grep 6:32561656

# PGEN
cat /exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/pgen/impute_recoded_selected_sample_filter_hq_var_6.pvar

# remove old results
rm -r results/meta/cojo/seq.13435.31/
rm results/meta/cojo/seq.13435.31_locus_chr6_31405620_32786291_ind_snps.tsv 
rm results/meta/logs/cojo/seq.13435.31.log

```

# Next-Flow pipeline for evaluation of the results
`/ssu/bsssu/dariush_coloc/testing_nf_pipeline`

- Permissions are granted to allow modifying the R scripts.

- Verifying the results using GLM model (Thu, 03:12, 18 to 20-Jun-23).


Dariush
