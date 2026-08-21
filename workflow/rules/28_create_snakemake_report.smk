rule create_snakemake_report:
    input:
        rules.QC_multiqc.output,
    output: 
        f"report_pipelineV{config['pipeline_version']}_{config['databases']['RVDB_description']}.html",
    log:
        "sn_logs/12_report/report.log",
    benchmark:
        "benchmarks/12_report/report_benchmark.txt"
    conda: "snakemake"
    params:
        snakefile=path_to_snakefile,
        sampledir=config["sampledir"],
    shell:
        """
        snakemake -s {params.snakefile} -c1 --config sampledir={params.sampledir} --report {output} \
        >{log} 2>&1
        """

