rule detect_illumina_chemistry:
    """
    Determine Illumina sequencer that was used to generate the data.
    This information (color-chemistry) will be used in trimming.
    """
    input:
        mate1="00_remove_contamination/{sample}_R1_001.tagged_filter.fastq.gz",
    output:
        chemistry="03_cutadapt/illumina_instrument_{sample}_R1_001.json"
    log:
        "sn_logs/03_cutadapt/illumina_instrument_{sample}.log",
    benchmark:
        "benchmarks/03_cutadapt/illumina_instrument_{sample}_benchmark.txt",
    threads: 4
    conda: "sm_python39"
    resources:
        mem_mb=200
    script:
          "../scripts/00_detect_illumina_chemistry/detect_illumina_chemistry.py"

rule trimming_paired:
    """
    Trim paired-end reads.
    """
    input:
        mate1="02_repair_mates/{sample}_R1_001.tagged_filter.repair.fastq",
        mate2="02_repair_mates/{sample}_R2_001.tagged_filter.repair.fastq",
        chemistry=rules.detect_illumina_chemistry.output.chemistry,
        illumina_adapters_R1=workflow.source_path(config["adapters"]["R1"]),
        illumina_adapters_R2=workflow.source_path(config["adapters"]["R2"]),
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
        quality_trim=lambda wildcards, input:
            get_quality_trim(input.chemistry),
    resources:
        mem_mb=200
    message: """--- Trimming on:
    {input}
    with cutadapt using {threads} cores and {resources.mem_mb} MB of RAM."""
    shell:
        """
        cutadapt \
        -a file:{input.illumina_adapters_R1} \
        -A file:{input.illumina_adapters_R2} \
        --times 1 \
        --cores {threads} \
        {params.quality_trim}=29 \
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
        chemistry=rules.detect_illumina_chemistry.output.chemistry,
        illumina_adapters_R1=workflow.source_path(config["adapters"]["R1"]),
    output:
        mate1_trimmed=temp("03_cutadapt/{sample}_R1_001.fastq"),
    log:
        "sn_logs/03_cutadapt/{sample}_trimmed_single.log",
    benchmark:
        "benchmarks/03_cutadapt/{sample}_trimmed_single_benchmark.txt",
    threads: 4
    conda: "sm_cutadapt"
    params:
        quality_trim=lambda wildcards, input:
            get_quality_trim(input.chemistry),
    resources:
        mem_mb=200
    message: """--- Trimming on:
    {input}
    using cutadapt on {threads} cores and {resources.mem_mb} MB of RAM."""
    shell:
        """
        cutadapt \
        -a file:{input.illumina_adapters_R1} \
        --times 1 \
        --cores {threads} \
        {params.quality_trim}=29 \
        --output {output.mate1_trimmed} \
        {input.mate1} > {log}
        """
