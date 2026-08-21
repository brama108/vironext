rule classify_contigs:
    """
    Classify de-novo assemblies.
    """
    input:
        rules.ML_predict_viruses.output.pred_fasta
    output:
        "13_ML_predict_virus/diamond/assembled_contigs_{sample}/assembled_contigs_{sample}_diamond.out"
    log:
        "sn_logs/13_ML_predict_virus/diamond/assembled_contigs_{sample}/assembled_contigs_{sample}_diamond.log"
    benchmark:
        "benchmarks/13_ML_predict_virus/diamond/assembled_contigs_{sample}/assembled_contigs_{sample}_diamond_benchmark.log"
    params:
        diamond_database=dynamic_path.databases.diamond__classify_contigs,
    conda:
        "sm_diamond"
    threads: 21
    resources:
        mem_mb=7680
    shell:
        """
        if [ -s {input} ]; then
            diamond blastx \
                --threads {threads} \
                -d {params.diamond_database} \
                --outfmt 6 qseqid sseqid stitle \
                --max-target-seqs 1 \
                -q {input} \
                --out {output}
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') - rule classify_contigs FAILED. Reason is probably empty or small input file" > {log}
            touch {output}
        fi
        """

