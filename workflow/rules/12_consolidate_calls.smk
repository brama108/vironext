rule consolidate_calls:
    """
    Consolidate primary and secondary classification. Only keep hits that are found in both bwa/micop and kaiju.
    """
    output:
         "14_validate_calls/final_classification/{sample}_R1R2_sorted_filtered_abundances_regions_micop_kaiju.txt"
    input:   
        micop_results=rules.extract_species_and_sort_abundances.output,
        kaiju_table=rules.summarize_results.output,
    benchmark:
        "benchmarks/14_validate_calls/final_classification/{sample}_R1R2_sorted_filtered_abundances_regions_micop_kaiju_benchmark.txt"
    params:
        dynamic_path.databases.taxonomy_db,
    threads: 1
    resources:
        mem_mb=5000
    conda: "sm_python39"
    script:
        "../scripts/04_consolidate_calls/consolidate_calls.py"