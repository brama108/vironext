rule extract_accessions:
    """
    Map taxids identified by micop to accessions in RVDB
    """
    output:
        "14_validate_calls/regions/{sample}_R1R2_sorted_filtered_abundances_regions.bed"
    input:    
        rules.extract_species_and_sort_abundances.output,
    benchmark:
        "benchmarks/14_validate_calls/regions/{sample}_R1R2_sorted_filtered_abundances_regions_benchmark.txt"
    params:
        dynamic_path.databases.taxonomy_db,
    threads: 1
    resources:
        mem_mb=600
    conda: "sm_python39"
    script:
        "../scripts/02_extract_accessions/extract_accessions.py"