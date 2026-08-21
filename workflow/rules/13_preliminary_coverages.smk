rule compute_preliminary_coverages:
    """
    Create coverage tables and graphs using samtools coverage
    """
    output:
        table="08_extract_high_coverage_species/{sample}_R1R2_sorted_filtered.txt",
    input:    
        rules.filter_by_mapQ.output,
    benchmark:
        "benchmarks/08_extract_high_coverage_species/{sample}_R1R2_sorted_filtered_benchmark.txt"
    conda:
        "sm_samtools",
    threads: 1
    resources:
        mem_mb=5120
    shell: """
        samtools coverage {input} > {output.table}
        """
         
