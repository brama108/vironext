rule validate_input:
    """
    Identifies empty fastq files.
    """
    output:
        "validate_input/{sample}_{readnum}_001_validate_input.out",
    input:
        lambda wildcards: f"{config["sampledir"]}/{wildcards.sample}_{wildcards.readnum}_001.fastq.gz",
    benchmark:
        "benchmarks/validate_input/{sample}_{readnum}_001_validate_input_benchmark.txt",
    threads: 1
    conda:  "sm_validatefastq"
    resources:
        mem_mb=1024
    shell:
        """
        biopet-validatefastq \
            --log_level warn \
            --fastq1 {input} \
            > {output} 2>&1
        # Check if output file exists AND has non-zero size
        if [ -s {output} ]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - rule validate_input_paired FAILED: please check {input}."
        fi
        """

rule QC_fastqc:
    """
    Perform quality control on untrimmed and trimmed
    """
    output: 
        html=multiext("01_fastqc/{when}/{sample}_{readnum}_001_{when}_fastqc", ".html", ".zip")
    input:
        fastq=rule_fastqc_get_input,
        input_validation=rules.validate_input.output,
    log:
        "sn_logs/01_fastqc/{when}/{sample}_{readnum}_001_fastqc.log"

    threads: 1
    conda:  "sm_fastqc"
    params:
        output_path="01_fastqc/{when}/"
    resources:
        mem_mb=300,
    message: """--- Quality check {wildcards.when} on:
        {input.fastq} 
    with Fastqc using {threads} cores and {resources.mem_mb} MB of RAM."""
    shell: """
        out_name=$(basename {input.fastq} | sed 's/\.fastq/_{wildcards.when}.fastq/')
        ln -fs {input.fastq} $out_name
        fastqc $out_name -o {params.output_path} \
        --dir {params} \
        2> {log}
        rm $out_name
        """

rule QC_multiqc:
    output:
        report("11_multiqc/multiqc_report.html",
            caption="../report/quality_control.txt", category="quality control"),
    input:
        fqc_html=expand("01_fastqc/{when}/{sample}_{readnum}_001_{when}_fastqc.html",
            sample=SAMPLE,
            readnum=READNUM,
            when=["before_trimming",
            "after_trimming"]),
        krona=expand("09_krona_plot/{sample}.html",
            sample=SAMPLE),
        wordclouds=expand("10_wordclouds/{sample}.html",
            sample=SAMPLE),
        table_html=expand("16_revalidation_blast/final_classification/{sample}.html",
            sample=SAMPLE),
        bigwigs=expand("05_coverage_tracks/bigwigs/{sample}_R1R2_coverage.bigwig",
            sample=SAMPLE),
        cov_graph=expand("15_realignment/05_coverages/{sample}.html",
            sample=SAMPLE),
        consensus=expand("15_realignment/07_consenus_sequence/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus.fasta",
            sample=SAMPLE),
        second_alignment_index=expand("15_realignment/04_mapping/sorted/{sample}_R1R2_sorted_filtered_selected_accessions_sorted.bam.bai",
            sample=SAMPLE),
        ML_report=expand("13_ML_predict_virus/report/{sample}_novel_virus_table.html",
            sample=SAMPLE),
        ML_contigs=expand("13_ML_predict_virus/report/{sample}_novel_virus_sequence.html",
            sample=SAMPLE),
        blast=expand(r"16_revalidation_blast/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_filtered.out",
            sample=SAMPLE),
        fai_index=expand("15_realignment/03_references/{sample}/{sample}_references.fasta.fai",
           sample=SAMPLE),
        # qm_metrics=expand("quality_metrics/{sample}.jpg",
        #    sample=SAMPLE),
        # qm_metrics_csv=expand("quality_metrics/{sample}.txt",
        #    sample=SAMPLE),
    threads: 1
    conda:
         "sm_multiqc"
    resources:
        mem_mb=250
    params:
        output_path="11_multiqc"
    shell: """
        multiqc --force --exclude kaiju --exclude htseq --outdir {params.output_path} .
        """
