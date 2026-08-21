rule generate_coverage_htmls:
    output:
        graph_html=report("15_realignment/05_coverages/{sample}.html",
            caption="../report/coverage_plots.txt",
            category="coverages"),
        # coverage_percent=temp("15_realignment/05_coverages/{sample}_coverage_percent.txt"),
        coverage_percent="15_realignment/05_coverages/{sample}_coverage_percent.txt",
    input:
        rules.generate_coverages.output,
        rules.generate_consensus.output,
        rules.create_final_identification_blast.output,
    benchmark:
        "benchmarks/15_realignment/05_coverages/{sample}_benchmark.txt"
    params:
        taxonomy_path=dynamic_path.databases.taxonomy_db,
        script_path_cov=script_path,
        script_path="../scripts/07_generate_htmls/generate_coverage_htmls.py",
    threads: 1
    resources:
        mem_mb=250
    conda: "sm_python39"
    script:
        "../scripts/08_generate_htmls/generate_coverage_htmls.py"
        # "{params.script_path}"
        ## "{script_path_cov}"

    #ToDo: uncomment the commented stuff and in common.smk delete the script_path variable. {params.script_path_cov} needd for rule generate_coverage_htmls. {script_path_cov} sowrks for report generation. params is needed for manual selection workflow
    # This is a snakemake bug (https://github.com/snakemake/snakemake/issues/1697)


rule generate_species_table_htmls:
    output:
        micop_filtered_table_html=report("16_revalidation_blast/final_classification/{sample}.html",
            caption="../report/micop_table.txt",
            category="identified taxon list"),
    input:
        rules.create_final_identification_blast.output,
        rules.generate_coverage_htmls.output.coverage_percent,
        rules.calculate_percentage_unambiguous_alignments.output,
    benchmark:
        "benchmarks/16_revalidation_blast/final_classification/{sample}_benchmark.txt"
    params:
        dynamic_path.databases.taxonomy_db,
        config["parameters"]["min_number_regions"],
        config["parameters"]["min_covbases"],
    threads: 1
    resources:
        mem_mb=250
    conda: "sm_python39"
    script:
        "../scripts/08_generate_htmls/generate_species_identification_table_htmls.py"