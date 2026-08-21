rule trimming_paired:
    """
    Trim paired-end reads.
    """
    input:
        mate1="02_repair_mates/{sample}_R1_001.tagged_filter.repair.fastq",
        mate2="02_repair_mates/{sample}_R2_001.tagged_filter.repair.fastq",
    output:
        mate1_trimmed=temp("03_cutadapt/{sample}_R1_001.fastq"),
        mate2_trimmed=temp("03_cutadapt/{sample}_R2_001.fastq"),
    log:
        "sn_logs/03_cutadapt/{sample}_trimmed_paired.log",
    benchmark:
        "benchmarks/03_cutadapt/{sample}_trimmed_paired_benchmark.txt",
    threads: 4
    conda: "sm_cutadapt"
    params:
        illumina_adapter_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        illumina_adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    resources:
        mem_mb=200
    message: """--- Trimming on:
    {input}
    with R1 illumina_adapter {params.illumina_adapter_R1} and R2 illumina_adapter {params.illumina_adapter_R2}
    with cutadapt using {threads} cores and {resources.mem_mb} MB of RAM."""
    shell:
        """
        cutadapt \
        -a {params.illumina_adapter_R1} \
        -A {params.illumina_adapter_R2} \
        --times 1 \
        --cores {threads} \
        --nextseq-trim=30 \
        --output {output.mate1_trimmed} \
        --paired-output {output.mate2_trimmed} \
        {input.mate1} {input.mate2} > {log}
        """

rule trimming_single:
    """
    Trim single-end reads.
    """
    input:
        mate1="00_remove_contamination/{sample}_R1_001.tagged_filter.fastq.gz",
    output:
        mate1_trimmed=temp("03_cutadapt/{sample}_R1_001.fastq"),
    log:
        "sn_logs/03_cutadapt/{sample}_trimmed_single.log",
    benchmark:
        "benchmarks/03_cutadapt/{sample}_trimmed_single_benchmark.txt",
    threads: 4
    conda: "sm_cutadapt"
    params:
        illumina_adapter_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    resources:
        mem_mb=200
    message: """--- Trimming on:
    {input}
    with R1 illumina_adapter {params.illumina_adapter_R1} 
    using cutadapt on {threads} cores and {resources.mem_mb} MB of RAM."""
    shell:
        """
        cutadapt \
        -a {params.illumina_adapter_R1} \
        --times 1 \
        --cores {threads} \
        --nextseq-trim=30 \
        --output {output.mate1_trimmed} \
        {input.mate1} > {log}
        """
