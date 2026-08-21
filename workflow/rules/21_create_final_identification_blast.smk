rule create_final_identification_blast:
    """"
    Use filtered blast output to create tje final species identification table.
    """
    input:
        rules.filter_identification_table.output,
        rules.evaluate_blast.output,
    output:
        "16_revalidation_blast/final_classification/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification.out",
    log:
        "sn_logs/16_revalidation_blast/final_classification/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification.log",
    benchmark:
        "benchmarks/16_revalidation_blast/final_classification/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification.log_benchmark.txt"
    params:
        blast_rvdb_database=dynamic_path.databases.blast_RVDB,
    conda: "sm_python39"
    script:
        "../scripts/10_create_final_identification_blast/create_final_identification_blast.py"