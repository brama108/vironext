use rule sort_and_convert as convert_sam with:
    output:
       temp("06_mapQ/{sample}_R1R2_sorted_filtered.bam"),
    input:    
        rules.filter_by_mapQ.output,
    benchmark:
        "benchmarks/06_mapQ/{sample}_R1R2_sorted_filtered_benchmark.txt"

use rule index_bam as index_first_alignment with:
    output:
       "06_mapQ/{sample}_R1R2_sorted_filtered.bam.bai",
    input:    
        rules.convert_sam.output,
    benchmark:
        "benchmarks/06_mapQ/{sample}_R1R2_sorted_filtered.bam_benchmark.txt"

rule extract_subalignment:
    """
    Extract regions from alignment file.
    """
    output:
         temp("14_validate_calls/extracted_reads/{sample}_R1R2_sorted_filtered_regions.bam"),
    input:    
        bed=rules.extract_accessions.output,
        alignment=rules.convert_sam.output,
        index_placeholder=rules.index_first_alignment.output,
    benchmark:
        "benchmarks/14_validate_calls/extracted_reads/{sample}_R1R2_sorted_filtered_regions_benchmark.txt"
    threads: 7
    resources:
        mem_mb=150
    conda: "sm_samtools"
    shell: """
        # Avoid crash if samtools returns with non-zero exit code (due to small file size). Error will be handled in the if statement
        set +e
        samtools view  {input.alignment} \
        -@ {threads} \
        --region-file {input.bed} -o {output}

        if [ $? -ne 0 ]; then
            echo "rule extract_subalignment FAILED. Reason is probably empty or small input file. Check logs"
            touch {output}
        fi
        """

rule extract_reads:
    """
    Extract fastq reads from alignment file.
    """
    output:
         temp("14_validate_calls/reads/{sample}_R1R2_sorted_filtered_regions.fastq")
    input:    
        rules.extract_subalignment.output,
    benchmark:
        "benchmarks/14_validate_calls/reads/{sample}_R1R2_sorted_filtered_regions_benchmark.txt"
    threads: 2
    resources:
        mem_mb=100
    conda: "sm_samtools"
    shell: """ 
        samtools fastq \
        -@ {threads} {input} > {output}
        """


        

        