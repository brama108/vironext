rule remove_contamination:
    """
    Remove host and other contaminating sequences
    """
    output:
        html_fastqscreen="00_remove_contamination/{sample}_{readnum}_001_screen.html",
        txt_fastqscreen="00_remove_contamination/{sample}_{readnum}_001_screen.txt",
        tagged="00_remove_contamination/{sample}_{readnum}_001.tagged.fastq.gz",
        filtered="00_remove_contamination/{sample}_{readnum}_001.tagged_filter.fastq.gz",

    input:
        input_validation=rules.validate_input.output,
        fastq=lambda wildcards: f"{config["sampledir"]}/{wildcards.sample}_{wildcards.readnum}_001.fastq.gz",
    benchmark:
        "benchmarks/00_remove_contamination/{sample}_{readnum}_001_benchmark.txt"
    params:
        output_path="00_remove_contamination",
        fastqscreen_config=workflow.source_path("../../config/fastq_screen.conf"),
    log:
        "sn_logs/00_remove_contamination/{sample}_{readnum}_001.log"
    threads: 6
    resources:
        mem_mb=12288
    conda:  "sm_fastq-screen"
    message: """--- Contamination check of raw data on:
    {input.fastq}
    with FastqScreen using {threads} cores and {resources.mem_mb} MB of RAM."""
    shell: """
        fastq_screen \
        --conf {params.fastqscreen_config} \
        --aligner bowtie2 \
        --threads {threads} \
        --outdir {params.output_path} \
        --nohits \
        --force \
        {input.fastq} \
        2> {log}
             """

rule repair_paired_ends:
    output: 
        mate1=temp("02_repair_mates/{sample}_R1_001.tagged_filter.repair.fastq"),
        mate2=temp("02_repair_mates/{sample}_R2_001.tagged_filter.repair.fastq"),
        singeltons="02_repair_mates/{sample}_001.tagged_filter.repair.singeltons.fastq",
    input:  
        mate1="00_remove_contamination/{sample}_R1_001.tagged_filter.fastq.gz",
        mate2="00_remove_contamination/{sample}_R2_001.tagged_filter.fastq.gz",
    log:
        "sn_logs/02_repair_mates/{sample}_R1R2_001_fastqscreen.log"
    threads: 2
    conda:  "sm_bbmap"
    resources:
        mem_mb=9216
    message: """--- Repairing mates:
    {input.mate1} and {input.mate2} 
    with repair.sh (bbtools) using {threads} cores and {resources.mem_mb} MB of RAM."""
    shell: """
        repair.sh \
        in={input.mate1} \
        in2={input.mate2} \
        -Xmx{resources.mem_mb}m \
        -out={output.mate1} \
        -out2={output.mate2} \
        -outs={output.singeltons} \
        2> {log}
        """
