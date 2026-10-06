rule extract_unmapped_reads_paired:
    """
    Extract reads that failed to align to viral genomes.
    """
    output:
        mate1_unmapped=temp("13_ML_predict_virus/{sample}_unmapped_R1.fastq"),
        mate2_unmapped=temp("13_ML_predict_virus/{sample}_unmapped_R2.fastq"),
    input:
        rules.sort_and_convert.output
    benchmark:
        "benchmarks/13_ML_predict_virus/extract_unmapped_reads_paired/{sample}_unmapped_benchmark.txt"
    conda:
        "sm_samtools"
    threads: 2
    resources:
        mem_mb=100
    shell: """
        samtools view {input} \
            | awk '$3 == "*" && length($10) > 50' \
            > 13_ML_predict_virus/{wildcards.sample}_unmapped.sam


        if [ -s 13_ML_predict_virus/{wildcards.sample}_unmapped.sam ]; then
            samtools fastq 13_ML_predict_virus/{wildcards.sample}_unmapped.sam \
                -1 {output.mate1_unmapped} \
                -2 {output.mate2_unmapped} \
                -0 /dev/null -s /dev/null -n
        else
            # if there are no unmapped reads then create empty files
            touch {output.mate1_unmapped}
            touch {output.mate2_unmapped}
        fi
        """

rule extract_unmapped_reads_single:
    """
    Extract reads that failed to align to viral genomes.
    """
    output:
        mate1_unmapped=temp("13_ML_predict_virus/{sample}_unmapped_R1.fastq"),
    input:
        rules.sort_and_convert.output
    benchmark:
        "benchmarks/13_ML_predict_virus/extract_unmapped_reads_single/{sample}_unmapped_benchmark.txt"
    conda:
        "sm_samtools"
    threads: 2
    resources:
        mem_mb=100
    shell: """
        samtools view {input} \
            | awk '$3 == "*" && length($10) > 50' \
            > 13_ML_predict_virus/{wildcards.sample}_unmapped.sam


        if [ -s 13_ML_predict_virus/{wildcards.sample}_unmapped.sam ]; then
            samtools fastq 13_ML_predict_virus/{wildcards.sample}_unmapped.sam > {output.mate1_unmapped}
        else
            # if there are no unmapped reads then create empty files
            touch {output.mate1_unmapped}
        fi
        """

rule assemble_contigs_paired:
    """
    Assemble contigs from reads that fail to align to viral genomes.
    """
    output:
        "13_ML_predict_virus/assembled_contigs_{sample}/contigs.fasta"
    input:
        mate1_unmapped=rules.extract_unmapped_reads_paired.output.mate1_unmapped,
        mate2_unmapped=rules.extract_unmapped_reads_paired.output.mate2_unmapped
    benchmark:
        "benchmarks/13_ML_predict_virus/assembled_contigs_{sample}/contigs_benchmark.txt"
    params:
        tmp_dir=config["tmp_dir"]
    conda:
        "sm_spades"
        #"../envs/spades.yml"
    threads: 10
    resources:
        mem_mb=17408
    shell: """
        if [ -s {input.mate1_unmapped} ]; then
            spades.py \
                --rnaviral \
                -t {threads} \
                --tmp-dir {params.tmp_dir} \
                -1 {input.mate1_unmapped} \
                -2 {input.mate2_unmapped} \
                -o 13_ML_predict_virus/assembled_contigs_{wildcards.sample}
        else
            touch {output}
        fi
        """

rule assemble_contigs_single:
    """
    Assemble contigs from reads that fail to align to viral genomes.
    """
    output:
        "13_ML_predict_virus/assembled_contigs_{sample}/contigs.fasta"
    input:
        mate1_unmapped=rules.extract_unmapped_reads_single.output.mate1_unmapped,
    benchmark:
        "benchmarks/13_ML_predict_virus/assembled_contigs_{sample}/contigs_benchmark.txt"
    params:
        tmp_dir=config["tmp_dir"]
    conda:
        "sm_spades"
        #"../envs/spades.yml"
    threads: 10
    resources:
        mem_mb=17408
    shell: """
        if [ -s {input.mate1_unmapped} ]; then
            spades.py \
                --rnaviral \
                -t {threads} \
                --tmp-dir {params.tmp_dir} \
                -s {input.mate1_unmapped} \
                -o 13_ML_predict_virus/assembled_contigs_{wildcards.sample}
        else
            touch {output}
        fi
        """

rule contig_mapping:
    """
    Mapping unclassified reads back to their contigs.
    """
    output:
        temp("13_ML_predict_virus/contigs_remap/{sample}_contig_remap.sam")
    input:
        contigs="13_ML_predict_virus/assembled_contigs_{sample}/contigs.fasta",
        reads=lambda wildcards: [
            f"13_ML_predict_virus/{wildcards.sample}_unmapped_R1.fastq"
        ] + ([f"13_ML_predict_virus/{wildcards.sample}_unmapped_R2.fastq"] if "R2" in READNUM else []),
    log:
        "sn_logs/13_ML_predict_virus/{sample}_contig_align.log"
    conda:
        "sm_bwa-mem2"
    threads: 5
    resources: mem_mb=5000
    shell: """
        # Build index from assembled contigs sequences
        mkdir -p 13_ML_predict_virus/contigs_bwa_index
        bwa-mem2 index -p 13_ML_predict_virus/contigs_bwa_index/{wildcards.sample} {input.contigs}

        # Map reads to coontig index
        bwa-mem2 mem -a -t {threads} 13_ML_predict_virus/contigs_bwa_index/{wildcards.sample} {input.reads} > {output} \
        2>{log}
        """

rule contig_mapping_stats:
    """
    Summarize remapping results
    """
    output:
        "13_ML_predict_virus/contigs_remap/remap_stats_{sample}.tsv"
    input:
        "13_ML_predict_virus/contigs_remap/{sample}_contig_remap.sam"
    conda:
        "sm_samtools"
    threads: 1
    resources: mem_mb=1000
    shell: """
        samtools sort {input} \
            | samtools coverage - \
            > {output}
        """


rule ML_predict_viruses:
    """
    Use VirHunter to predict if contigs originated from a virus.
    """
    input:
         "13_ML_predict_virus/assembled_contigs_{sample}/contigs.fasta",
    output:
        pred_csv="13_ML_predict_virus/predicted_viruses_{sample}/contigs_predicted.csv",
        pred_fasta="13_ML_predict_virus/predicted_viruses_{sample}/contigs_viral.fasta"
    benchmark:
        "benchmarks/13_ML_predict_virus/predicted_viruses_{sample}/contigs_predicted_benchmark.txt"
    params:
        vh_path=workflow.source_path("../scripts/virhunter/predict.py"),
        vh_utils=workflow.source_path("../scripts/virhunter/utils/preprocess.py"),
        vh_mod10=workflow.source_path("../scripts/virhunter/models/model_10.py"),
        vh_mod7=workflow.source_path("../scripts/virhunter/models/model_7.py"),
        vh_mod5=workflow.source_path("../scripts/virhunter/models/model_5.py"),
        vh_model=dynamic_path.databases.virhunter_weights,
        vh_config=workflow.source_path("../../config/config_virhunter.yml")
    conda:
        "sm_virhunter"
        #"../envs/virhunter.yml"
    threads: 25
    resources:
        #partition="GPU"
        mem_mb=250000
    shell: """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/

        sed -e 's|<G>|{input}|' \
            -e 's|<H>|{params.vh_model}|' \
            -e 's|<I>|13_ML_predict_virus/predicted_viruses_{wildcards.sample}|' \
            {params.vh_config} > 13_ML_predict_virus/config_{wildcards.sample}.yml

        if [ -s {input} ]; then
            python {params.vh_path} 13_ML_predict_virus/config_{wildcards.sample}.yml
        else
            touch {output.pred_csv}
            touch {output.pred_fasta}
        fi
        """
        
