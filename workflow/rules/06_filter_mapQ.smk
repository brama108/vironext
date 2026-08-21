rule filter_by_mapQ:
    output:
        temp("06_mapQ/{sample}_R1R2_sorted_filtered.sam")
    input:
        sorted=rules.sort_and_convert.output,
        index=rules.index_bam.output
    log:
        "sn_logs/06_mapQ/{sample}_R1R2_001_filter_by_mapQ.log"
    benchmark:
        "benchmarks/06_mapQ/{sample}_R1R2_001_filter_by_mapQ_benchmark.txt"
    conda:
        "sm_sambamba",
    threads:11
    resources:
        mem_mb=1024
    shell: """
        sambamba view -t {threads} -h -f sam \
        -F "(mapping_quality<1 or mapping_quality>=20)" \
        {input.sorted} \
        -o {output} \
        >{log} 2>&1
       """