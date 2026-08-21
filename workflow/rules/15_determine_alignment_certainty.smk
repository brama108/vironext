rule get_number_unambiguous_alignments:
    """
    Alignments with mapQ >= 30 are considered unambiguous. Alignment counts are calculated per selected accession.
    """
    output:
        f"{path_to_abiguous_alignments}/{{sample}}/{{accession}}_unambig.txt",
    input:    
        sorted=rules.process_alignment.output,
        index=rules.index_second_alignment.output,
    log:
        f"sn_logs/{path_to_abiguous_alignments}/{{sample}}/{{sample}}_{{accession}}_unambiguous_alignments.txt.log"
    benchmark:
         f"benchmarks/{path_to_abiguous_alignments}/{{sample}}/{{sample}}_{{accession}}_unambiguous_alignments_benchmark.txt",
    conda:
        "sm_sambamba",
    threads:11
    resources:
        mem_mb=1024
    shell: """
        sambamba view -t {threads} -f sam \
        -F "(mapping_quality>=60)" \
        {input.sorted} \
        {wildcards.accession} \
        | wc -l > {output}
        >{log} 2>&1
       """

rule get_total_number_alignments:
    """
    Count total alignments per selected accession.
    """
    input:
        sorted=rules.process_alignment.output,
        index=rules.index_second_alignment.output,
    output:
        f"{path_to_abiguous_alignments}/{{sample}}/{{accession}}_total_number_alignments.txt"
    log:
        f"sn_logs/{path_to_abiguous_alignments}/{{sample}}/{{accession}}_total_number_alignments.log"
    benchmark:
         f"benchmarks/{path_to_abiguous_alignments}/{{sample}}/{{sample}}_{{accession}}_total_number_alignments.txt",
    threads: 11
    resources:
        mem_mb=1024
    conda: "sm_sambamba"
    shell:
        """
        sambamba view -t {threads} -f sam \
        {input.sorted} \
        {wildcards.accession} \
        | wc -l > {output}
        >{log} 2>&1
        """

rule calculate_percentage_unambiguous_alignments:
    input:
        total=get_total_files,
        unambig=get_unambig_files,
    output:
        f"{path_to_abiguous_alignments}/{{sample}}/{{sample}}_pct.txt",
    params: accession_dir="15_realignment/02_max_reads_accessions/{sample}/accession_files",
    log:
        f"sn_logs/{path_to_abiguous_alignments}/{{sample}}/{{sample}}_pct.log"
    benchmark:
        f"benchmarks/{path_to_abiguous_alignments}/{{sample}}/{{sample}}_pct_benchmark.txt"
    conda: "sm_python39"
    threads: 1
    script: "../scripts/06_calculate_percentage_unamgiguous_alignments/calculate_percentage_unambiguous_alignments.py"