rule mapping:
    """
    Mapping cleaned reads. Depending on read type either SE or PE data serves as input.
    """
    output:
        temp("04_bwa_micop/{sample}_R1R2.sam"),
    input:
        lambda wildcards: [
            f"03_cutadapt/{wildcards.sample}_R1_001.fastq"
        ] + ([f"03_cutadapt/{wildcards.sample}_R2_001.fastq"] if "R2" in READNUM else []),
    params:
        virus_database=dynamic_path.databases.condensed_RVDB,
        read_group="{sample}",
    log:
        "sn_logs/04_bwa_micop/{sample}_R1R2_001_bwa.log"
    benchmark:
        "benchmarks/04_bwa_micop/{sample}_R1R2_001_bwa_benchmark.txt"
    conda:
        "sm_bwa-mem2"
    threads: 26
    resources: mem_mb=43008
    shell: """
        bwa-mem2 mem -a -t {threads} -R "@RG\\tID:{wildcards.sample}\\tSM:{params.read_group}" {params.virus_database} {input} > {output} \
        2>{log}
        """