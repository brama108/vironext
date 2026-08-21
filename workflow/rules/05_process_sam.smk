rule sort_and_convert:
    """
    Sort sam and convert to bam.
    """
    output:
        temp("04_bwa_micop/sorted/{sample}_R1R2_sorted.bam"),
    input:
        rules.mapping.output
    benchmark:
        "benchmarks/04_bwa_micop/sorted/{sample}_R1R2_sorted_benchmark.txt"
    threads: 8
    resources:
        mem_mb=8192
    conda: "sm_samtools"
    shell: """
        samtools sort -@ {threads} \
        -O bam \
        -o {output} \
        {input}
        """

rule index_bam:
    output:
        temp("04_bwa_micop/sorted/{sample}_R1R2_sorted.bam.bai"),
    input:    
        rules.sort_and_convert.output
    conda:
        "sm_samtools",
    threads: 2
    resources:
        mem_mb=100
    shell: """
        samtools index -@{threads} {input}
        """