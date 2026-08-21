rule create_novel_virus_report:
    """
    Generate html table with virus predition
    """
    input:
        rules.ML_predict_viruses.output.pred_csv,
        rules.classify_contigs.output,
        rules.contig_mapping_stats.output
    output:
        report("13_ML_predict_virus/report/{sample}_novel_virus_table.html",
            category="novel virus prediction",
            subcategory="results table")
    benchmark:
        "benchmarks/13_ML_predict_virus/report/{sample}_novel_virus_table_benchmark.txt"
    params:
        readcutoff=dynamic_path.parameters.read_cutoff_ML,
    conda:
        "sm_python39"
    threads: 1
    resources:
        mem_mb=500
    script:
        "../scripts/08_generate_htmls/generate_novel_virus_table.py"

rule create_virus_contig_report:
    """
    Generate html version of .fasta contigs
    """
    input:
        rules.ML_predict_viruses.output.pred_fasta
    output:
        report("13_ML_predict_virus/report/{sample}_novel_virus_sequence.html",
            category="novel virus prediction",
            subcategory="sequences")
    benchmark:
        "benchmarks/13_ML_predict_virus/report/{sample}_novel_virus_sequence_benchmark.txt"
    conda:
        "sm_python39"
    threads: 1
    resources:
        mem_mb=500
    shell:
        """
        echo '<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
  </head>
  <body>
      <pre style="line-height: 1; font-family: monospace, monospace !important; font-size: 1em; white-space: pre; color: DarkBlue">
      ' > {output}

        echo "<h1>Sample: {wildcards.sample}</h1>" >> {output}
        echo "$(< {input})" >> {output}
        echo "
      </pre>
  </body>
</html>" >> {output}
        """