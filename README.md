# pqtl_pipeline_finemap
Fine mapping analysis within the pQTL pipeline project


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
ls -1 /exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats_digits_not_flipped/output/seq*/seq*.gwaslab.tsv.bgz > conf/path_meta_all.txt
```

# mapping files
`/exchange/healthds/pQTL/results/INTERVAL/qced_sumstats/table.snp_mapping.tsv.gz`

- Verifying the results using GLM model (Thu, 03:12, 18 to 20-Jun-23).

- Lunching the locus breaker and COJO on the latest meta-analysis results (Mon, 23:55, 24-Jun-24)

- Running locus definition and COJO using updated PGEN/BED genotypes got finished (Tue, 19:30, 10-Jul-24).   

Dariush
