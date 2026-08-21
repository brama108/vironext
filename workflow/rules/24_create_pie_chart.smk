rule create_pie_chart:
    """
    Create pie chart depicting the relative abundance of identified viruses
    """
    output:
        report("09_krona_plot/{sample}.html", caption="../report/pie_chart_caption.txt", category="relative abundance", subcategory="pie chart"),
    input:
        rules.create_final_identification_blast.output,
    log:
        "sn_logs/09_krona_plot/{sample}_R1R2_001_generate_krona_plot.log"
    benchmark:
        "benchmarks/09_krona_plot/{sample}_R1R2_001_generate_krona_plot_benchmark.txt"
    params:
       krona_taxonomy=dynamic_path.databases.krona_taxonomy,
    conda:
        "sm_krona"
    resources:
        mem_mb=250
    shell: """
        ktImportTaxonomy \
        -t 1 \
        -m 3 \
        -o {output} \
        -tax {params.krona_taxonomy} \
        <(awk -F"\t" '$8 == "True"' {input}) \
        >{log}
        """
