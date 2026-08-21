rule generate_coverages:
    """
    Create coverage graphs using samtools coverage
    """
    output:
        "15_realignment/05_coverages/{sample}_R1R2_sorted_filtered_selected_accessions_sorted.graph.txt",
    input:    
        rules.process_alignment.output,
    benchmark:
        "benchmarks/15_realignment/05_coverages/{sample}_R1R2_sorted_filtered_selected_accessions_sorted.graph_benchmark.txt"
    conda:
        "sm_samtools",
    threads: 1
    resources:
        mem_mb=120
    shell: """
        samtools coverage \
        --histogram \
        --ascii \
        --excl-flags 1536 \
        {input} \
        > {output}
        """
