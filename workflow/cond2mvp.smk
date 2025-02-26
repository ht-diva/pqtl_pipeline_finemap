# Summary of Snakemake Liftover Rule Development

# 1. **Initial Rule Setup:**
#    - Created a Snakemake rule that processes all VCF files from a specified directory.
# 2. **Dynamic File Discovery:**
#    - Implemented a function `get_vcf_files` to dynamically read all VCF files from the input directory.
# 3. **Wildcard Management:**
#    - Moved the renaming of VCF files (removal of `.vcf` extension) to the `get_vcf_files` function for cleaner code.
# 4. **Output Format Change:**
#    - Updated the rule to convert VCF files to TXT format after the liftover process.
# 5. **Final Snakemake Rules:**
#    - `rule all` ensures that all lifted TXT files are generated.
#    - `rule liftover_vcf` runs the liftover tool using a shell command and saves the output as TXT.



# read the configuration file
configfile: "conf/config_mvp.yml"

from pathlib import Path
import os

# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("workspace"), file_path))

def get_vcf_files(input_dir):
   """Read all VCF files from the specified directory."""
   """and remove the .vcf extension."""
   return [f.replace('.vcf', '') for f in os.listdir(input_dir) if f.endswith(".vcf")]


# Discover VCF files in the input directory and save lis of them
vcf_directory = ws_path("VCF")
vcf_files = get_vcf_files(vcf_directory)


rule all:
    input:
        ws_path("mr_results_unique_associations.tsv"),
        expand(ws_path("VCF_lifted/{filename}.txt"), filename=vcf_files),
        expand(ws_path("mvp/{filename}.sentinel"), filename=vcf_files)


# Remove multiple associations for variants in MR results and 
# for each variant, extract positions of variants in the locus from conditional data
# then, save positions and alleles in VCF format
rule grep_unique_targets:
    input:
        config.get("path_mvp")
    output:
        unique = ws_path("mr_results_unique_associations.tsv"),
        annotated = ws_path("mr_results_annotated.tsv")
    params:
        config.get("workspace")
    log:
        ws_path("logs/mr_results_unique_associations.log")
    conda:
        "envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    script:
        "scripts/tidy_mr_results.R"


# Snakemake rule to lift over VCF files from GRCh37 to GRCh38
rule liftover_vcf:
    input:
        ws_path("VCF/{filename}.vcf"),
        ws_path("mr_results_unique_associations.tsv")
    output:
        ws_path("VCF_lifted/{filename}.txt")
    log:
        ws_path("logs/liftover/{filename}.log"),
    params:
        chain_file = config.get("chain_file"),
        hg37 = config.get("fasta_hg37"),
        hg38 = config.get("fasta_hg38"),
        script="workflow/scripts/liftover_bcftool.sh"
    shell:
        """
        {params.script} {input[0]} {output} {params.chain_file} {params.hg37} {params.hg38} 
        """

# rule to extract conditional data from RDS and merge them with lifted positions
rule extract_conditional:
    input:
        mvp = ws_path("mr_results_unique_associations.tsv"),
        pos38 = ws_path("VCF_lifted/{filename}.txt")
    output:
        sentinel=touch(ws_path("mvp/{filename}.sentinel"))
    params:
        path = ws_path("mvp"),
        script="workflow/scripts/extract_conditional.R"
    conda:
        "envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 60 + attempt * 30
    shell:
        """
        Rscript {params.script} {input.mvp} {input.pos38} {params.path}
        """
