rule run_validation:
    """
    Run kaiju.
    """
    output:
         "14_validate_calls/kaiju/{sample}_R1R2_sorted_filtered_regions_kaiju.txt"
    input:    
        rules.extract_reads.output,
    benchmark:
        "benchmarks/14_validate_calls/kaiju/{sample}_R1R2_sorted_filtered_regions_kaiju_benchmark.txt"
    params:
        kaiju_nodes_file = dynamic_path.databases.kaiju_nodes_dmp,
        kaiju_db_index_file = dynamic_path.databases.kaiju_index,
    threads: 28
    resources:
        mem_mb=102400
    conda: "sm_kaiju"
    shell: """
        kaiju -t {params.kaiju_nodes_file} \
            -f {params.kaiju_db_index_file} \
            -i {input} \
            -z {threads} \
            -o {output}
        """
 

rule summarize_results:
    """
    Summarize kaiju results.
    """
    output:
         "14_validate_calls/kaiju/{sample}_R1R2_sorted_filtered_regions_kaiju_table.tsv"
    input:    
        rules.run_validation.output,
    benchmark:
        "benchmarks/14_validate_calls/kaiju/{sample}_R1R2_sorted_filtered_regions_kaiju_table_benchmark.txt"
    params:
        kaiju_nodes_file = dynamic_path.databases.kaiju_nodes_dmp,
        kaiju_names_file = dynamic_path.databases.kaiju_names_dmp,
    threads: 4
    resources:
        mem_mb=5000
    conda: "sm_kaiju"
    shell: """
        kaiju2table \
            -t {params.kaiju_nodes_file} \
            -n {params.kaiju_names_file} \
            -r species  \
            -e \
            -o {output} \
            {input} \
        """
