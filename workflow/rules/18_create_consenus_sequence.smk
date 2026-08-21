rule generate_consensus:
    """
    Generates consensus sequence for each identified species
    """
    output:
        "15_realignment/07_consenus_sequence/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus.fasta",
    input:    
        rules.process_alignment.output,
    benchmark:
        "benchmarks/15_realignment/07_consenus_sequence/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_benchmark.txt"
    conda:
        "sm_samtools",
    threads: 1
    resources:
        mem_mb=2048
    shell: """
        samtools consensus --ambig --min-depth 5 --format fasta --output {output} {input}
        """

        #TD: update samtools to >=1.16 for --min-BQ 20