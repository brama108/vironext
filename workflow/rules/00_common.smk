import dynamic_yaml
import os
import sys

##### Helper functions #####
def load_config(config_path):
    "Use dynamic_yaml to resolve config paths."
    with open(config_path) as fileobj:
        cfg = dynamic_yaml.load(fileobj)
    return cfg
    
def rule_fastqc_get_input(wildcards):
    if (wildcards.when == "before_trimming"):
        input_file = f"{path}/{{sample}}_{{readnum}}_001.fastq.gz"
        trim = ""
    elif (wildcards.when == "after_trimming"):
        input_file = "03_cutadapt/{sample}_{readnum}_001.fastq"
        trim = "trimmed_"
    return input_file

# input function for rules get_number_unambiguous_alignments and get_total_number_alignments, return paths to all files produced by the checkpoint select_accessions
def get_total_files(wildcards):
    checkpoint_out_dir = checkpoints.select_accessions.get(sample=wildcards.sample).output.accession_files
    accession_ids = glob_wildcards(os.path.join(checkpoint_out_dir, "{accession}.txt")).accession
    return expand(f"{path_to_abiguous_alignments}/{{sample}}/{{accession}}_total_number_alignments.txt",
                  sample=wildcards.sample, accession=accession_ids)

def get_unambig_files(wildcards):
    checkpoint_out_dir = checkpoints.select_accessions.get(sample=wildcards.sample).output.accession_files
    accession_ids = glob_wildcards(os.path.join(checkpoint_out_dir, "{accession}.txt")).accession
    return expand(f"{path_to_abiguous_alignments}/{{sample}}/{{accession}}_unambig.txt",
                  sample=wildcards.sample, accession=accession_ids)

def extract_path_to_snakefile():
    if "-s" in sys.argv:
        index = sys.argv.index("-s")
    elif "--snakefile" in sys.argv:
        index = sys.argv.index("--snakefile")
    snakefile_path = sys.argv[index + 1]
    return snakefile_path

path = config["sampledir"]
wildcard_pattern = f"{path}/{{sample}}_{{readnum}}_001.fastq.gz"
SAMPLE, READNUM = glob_wildcards(wildcard_pattern)
read_type = "single_end" if "R2" not in READNUM else "paired_end"
path_to_abiguous_alignments = "15_realignment/06_ambiguous_alignment_stats"

report: workflow.source_path("../report/workflow.rst"),
SAMPLE = set(SAMPLE)
READNUM = set(READNUM)
path_to_snakefile = extract_path_to_snakefile()
script_path=("/").join((path_to_snakefile.rsplit("/",1)[0],"scripts/07_generate_htmls/generate_coverage_htmls.py")),
cfg_path = os.path.abspath(workflow.source_path('../../config/config.yaml')).split("file")[1]
dynamic_path = load_config(cfg_path)
