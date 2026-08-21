rule calculate_abundances:
    """
    Calculate abundances using MiCoP's probabilistic classification approach
    """
    output:
        "07_calculate_abundances/{sample}_R1R2_sorted_filtered_abundances.txt"
    input:
        rules.filter_by_mapQ.output
    params:
        compute_abundances_MiCoP=workflow.source_path("../scripts/01_compute_abundances/compute-abundances.py"),
        micop_acc2info=dynamic_path.databases.micop_acc2info,
    log:
        "sn_logs/07_calculate_abundances/{sample}_R1R2_001_calculate_abundances.log"
    benchmark:
        "benchmarks/07_calculate_abundances/{sample}_R1R2_001_calculate_abundances_benchmark.txt"
    conda:
        "sm_python39",
    threads: 1
    resources:
        mem_mb=250000
    shell: """
        python {params.compute_abundances_MiCoP} {input} \
        --virus \
        --read_cutoff 0 \
        --pct_id 0.6 \
        --acc2info {params.micop_acc2info} \
        --raw_counts \
        --output {output} \
        >{log}
        """

rule extract_species_and_sort_abundances:
    """
    Extract species level classifications from the MiCoP output file
    """
    output:
        "07_calculate_abundances/{sample}_R1R2_sorted_filtered_abundances_species.txt",
    input:
        rules.calculate_abundances.output,    
    benchmark:
        "benchmarks/07_calculate_abundances/{sample}_R1R2_sorted_filtered_abundances_species_benchmark.txt"
    conda:
        "sm_python39",
    threads: 1
    resources:
        mem_mb=250
    shell: """
        if [[ $(wc -l < {input}) -lt 2 ]]; then
            touch {output}
        else
            grep species {input} | sort -t $'\\t' -grk5 | cut -f1,4- | cut -d "|" -f1,7- | sed 's/Viruses|//g'> {output}
        fi
        """
