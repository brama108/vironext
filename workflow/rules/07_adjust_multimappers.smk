#TD: -f bam. späteres umwandeln in bam sortieren und indizieren unnötig.
rule adjust_multimappers:
    """
    Filter for uniquely mapping reads. This rule ensures that identified species
    with multimapping reads only will not be considered.
    This rule effectively reduces the number of FPs species passing filters.
    """
    input:
        rules.filter_by_mapQ.output,
    output:
        temp("06_mapQ/02_unique_mapper/{sample}_R1R2_sorted_filtered_unique.sam"),
    log:
        "sn_logs/02_unique_mapper/{sample}_R1R2_sorted_filtered_unique.txt",
    benchmark:
        "benchmarks/06_mapQ/{sample}_R1R2_001_filter_by_mapQ_benchmark.txt"
    conda:
        "sm_sambamba",
    threads: 2
    resources:
        mem_mb=150
    shell: """
        sambamba view -t {threads} -h --sam-input -f sam \
        -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
        {input} \
        -o {output} \
        >{log} 2>&1
       """