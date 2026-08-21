rule revalidate_blast:
    """"
    Blast consensus sequences against C-RVDB to add another layer of validation.
    """
    input:
        rules.generate_consensus.output,
    output:
        "16_revalidation_blast/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast.out",
    log:
        "sn_logs/16_revalidation_blast/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast.log"
    benchmark:
        "benchmarks/16_revalidation_blast/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_benchmark.txt"
    params:
        blast_rvdb_database=dynamic_path.databases.blast_RVDB,
    conda: "sm_blastn"
    threads: 4
    resources:
        mem_mb=4096
    shell: """
        blastn \
        -num_threads {threads} \
        -db {params.blast_rvdb_database} \
        -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -query {input} \
        -out {output} """

rule evaluate_blast:
    """
    Filter out sequences with conflicting species identifications.
    """
    input:
        rules.revalidate_blast.output,
    output:
        "16_revalidation_blast/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_filtered.out",
    benchmark:
        "benchmarks/16_revalidation_blast/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_filtered_benchmark.txt"
    conda: "sm_python39"
    params:
        taxonomy_path=dynamic_path.databases.taxonomy_db,
    script:
        "../scripts/09_evaluate_blast/evaluate_blast.py"
