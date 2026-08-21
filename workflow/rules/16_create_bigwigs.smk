rule extract_genome_sizes:
    output:
        temp("05_coverage_tracks/genome_sizes/{sample}_R1R2_sorted_filtered_selected_accessions_genome_sizes.txt"),
    input:
        rules.realign_reads.output,
    log:
        "sn_logs/05_coverage_tracks/genome_sizes/{sample}_R1R2_sorted_filtered_selected_accessions_genome_sizes.log",
    benchmark:
        "benchmarks/05_coverage_tracks/genome_sizes/{sample}_R1R2_sorted_filtered_selected_accessions_genome_sizes_references_benchmark.txt"
    threads: 1
    conda:
        "sm_python39",
    resources:
        mem_mb=100
    shell: """
        awk -v OFS="\t" '$1 ~ /^@SQ/ {{split($2, chr, ":") 
                                  split($3, size, ":") 
                                  print chr[2], size[2]}}' {input} > {output}
        """

rule create_bedgraph:
    output:
        temp("05_coverage_tracks/bedgraphs/{sample}_R1R2_sorted_filtered_selected_accessions_coverage.bdg"),
    input:
        rules.process_alignment.output,
    log:
        "sn_logs/05_coverage_tracks/bedgraphs/{sample}_R1R2_sorted_filtered_selected_accessions_create_bedgraph.log",
    benchmark:
        "benchmarks/05_coverage_tracks/bedgraphs/{sample}_R1R2_sorted_filtered_selected_accessions_create_bedgraph_benchmark.txt"
    conda:
        "sm_bedtools"
    threads: 1
    resources:
        mem_mb=100
    shell: """
    # Avoid crash if bedtools returns with non-zero exit code (due to small file size). Error will be handled in the if statement
    set +e
    mkdir -p 05_coverage_tracks/bedgraphs/{wildcards.sample} # create directory for storing temporary files
    bedtools genomecov -ibam {input} \
                       -bg \
        | LC_COLLATE=C sort -k1,1 -k2,2n \
        -T 05_coverage_tracks/bedgraphs/{wildcards.sample}/ \
        > {output}
    rm -r 05_coverage_tracks/bedgraphs/{wildcards.sample} # remove directory for storing temporary files
    if [ $? -ne 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - rule create_bedgraph FAILED. Reason is probably empty or small input file. Check logs" > {log}
        touch {output}
    fi 
    """

rule bedgraphtobigwig:
    output:
        "05_coverage_tracks/bigwigs/{sample}_R1R2_coverage.bigwig",
    input:
        bedgraph=rules.create_bedgraph.output,
        chrom_sizes=rules.extract_genome_sizes.output,
    log:
        "sn_logs/05_coverage_tracks/bigwigs/{sample}_R1R2_bedGraphToBigWig.log",
    benchmark:
        "benchmarks/05_coverage_tracks/bigwigs/{sample}_R1R2_bedGraphToBigWig_benchmark.txt"
    conda:
        "sm_bedgraphtobigwig"
    threads: 1
    resources:
        mem_mb=100
    shell: """
    if [ -s {input.bedgraph} ]; then
        bedGraphToBigWig \
            {input.bedgraph} \
            {input.chrom_sizes} \
            {output}
            >{log}
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') - rule bedgraphtobigwig FAILED. Reason is probably empty or small input file" > {log}
        touch {output}
    fi
    """
