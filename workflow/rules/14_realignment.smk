# ToDo: The rules of this file could go into a module/subworkflow
rule alignment2regions:
    """
    Convert bam to bed file.
    """
    output:
         "15_realignment/01_regions/{sample}_R1R2_sorted_filtered.bed",
    input:    
        bam=rules.adjust_multimappers.output,
    benchmark:
        "benchmarks/15_realignment/01_regions/{sample}_R1R2_sorted_filtered_benchmark.txt"
    threads: 1
    resources:
        mem_mb=100
    conda: "sm_bedtools"
    shell: """
        bamToBed -i {input} > {output}
        """

rule merge_regions:
    """
    Merge overlapping regions.
    """
    output:
        "15_realignment/01_regions/{sample}_R1R2_sorted_filtered_merged.bed",
    input:
        rules.alignment2regions.output,
    benchmark:
        "benchmarks/15_realignment/01_regions/{sample}_R1R2_sorted_filtered_merged_benchmark.txt"
    threads: 1
    resources:
        mem_mb=51200
    conda: "sm_bedtools"
    shell: """ 
        bedtools merge -i {input} > {output}
        """

checkpoint select_accessions:
    """
    Find accessions with maximum number of reads per species.
    Apply filter for non-overlapping genome regions and covered bases.
    """
    output:
        selected_accessions="15_realignment/02_max_reads_accessions/{sample}/{sample}_R1R2_sorted_filtered_selected_accessions.txt",
        selected_species_taxids="15_realignment/02_max_reads_accessions/{sample}/{sample}_R1R2_sorted_filtered_selected_species_taxids.txt",
        accession_files=directory("15_realignment/02_max_reads_accessions/{sample}/accession_files"),
    input:
        merged_regions=rules.merge_regions.output,
        consolidated_species_list=rules.extract_species_and_sort_abundances.output,
        coverage_table=rules.compute_preliminary_coverages.output.table,
    benchmark:
        "benchmarks/15_realignment/02_max_reads_accessions/{sample}_R1R2_sorted_filtered_selected_merged_benchmark.txt"
    params:
        dynamic_path.databases.taxonomy_db,
        config["parameters"]["min_number_regions"],
        config["parameters"]["min_covbases"],
        config["parameters"]["min_region_length"],
    threads: 1
    resources:
        mem_mb=50
    conda: "sm_python39"
    script:
        "../scripts/03_select_accessions/select_accessions.py"

rule filter_identification_table:
    """
    Generates an identification table based on secondary classification and coverage filtering.
    """
    output:
        "14_validate_calls/filtered_identifications/{sample}_R1R2_sorted_filtered_abundances_regions_micop_kaiju_filtered.txt"
    input:
        path_to_micop_output = rules.extract_species_and_sort_abundances.output,
        micop_kaiju_output=rules.consolidate_calls.output,
        selected_species_taxids_path = rules.select_accessions.output.selected_species_taxids,
        #filtered_blast_output = rules.evaluate_blast.output,
    benchmark:
        "benchmarks/14_validate_calls/filtered_identifications/{sample}_R1R2_sorted_filtered_abundances_regions_micop_kaiju_filtered_benchmark.txt"
    params:
        dynamic_path.databases.taxonomy_db,
    threads: 1
    resources:
        mem_mb=5000
    conda: "sm_python39"
    script:
        "../scripts/05_filter_identification_table/filter_identification_table.py"

rule extract_references:
    """
    Extract sequences of the selected accessions and write them to a multifasta file.
    """
    output:
        "15_realignment/03_references/{sample}/{sample}_references.fasta",
    input:
        rules.select_accessions.output.selected_accessions,
    benchmark:
        "benchmarks/15_realignment/03_references/{sample}/{sample}_references_benchmark.txt"
    params:
        cdbtools_index=dynamic_path.databases.cdb_index,
    threads: 1
    resources:
        mem_mb=5000
    conda: "sm_cdbtools"
    shell: """
        cat {input} |
        cdbyank  {params.cdbtools_index}  \
        > {output}
        """

rule create_mapping_index:
    """
    Index multifasta reference file for bwa.
    """
    output:
        multiext("15_realignment/03_references/{sample}/{sample}_references.fasta", ".ann", ".amb", ".0123", ".pac", ".bwt.2bit.64"),
    input:
        rules.extract_references.output,
    benchmark:
        "benchmarks/15_realignment/03_references/{sample}/{sample}_references_bam_index_benchmark.txt"
    threads: 1
    resources:
        mem_mb=5000
    conda: "sm_bwa-mem2"
    shell: """
    bwa-mem2 index {input}
    """

rule create_fasta_index:
    """
    Index multifasta reference file or IGV.
    """
    output:
       "15_realignment/03_references/{sample}/{sample}_references.fasta.fai",
    input:
        rules.extract_references.output,
    benchmark:
        "benchmarks/15_realignment/03_references/{sample}/{sample}_references_fasta_index_benchmark.txt"
    threads: 1
    resources:
        mem_mb=5000
    conda: "sm_samtools"
    shell: """
    if [ -s {input} ]; then
        samtools faidx {input}
    else
        touch {output}
    fi
    """

rule realign_reads:
    """
    Map trimmed reads to selected sequences. Input depends on whether SE or PE reads are used.
    """
    output:
        temp("15_realignment/04_mapping/{sample}_R1R2_sorted_filtered_selected_accessions.sam"),
    input:
        mates=lambda wildcards: [
            f"03_cutadapt/{wildcards.sample}_R1_001.fastq"
        ] + ([f"03_cutadapt/{wildcards.sample}_R2_001.fastq"] if "R2" in READNUM else []),
        index=rules.create_mapping_index.output,
        reference=rules.extract_references.output,
    benchmark:
        "benchmarks/15_realignment/04_mapping/{sample}_R1R2_sorted_filtered_selected_accessions_benchmark.txt"
    params:
        read_group="{sample}",
    threads: 16
    resources:
        mem_mb=12288
    conda: "sm_bwa-mem2"
    shell: """
        bwa-mem2 mem -t {threads} \
        {input.reference} \
        -R "@RG\\tID:{wildcards.sample}\\tSM:{params.read_group}" \
        {input.mates} \
        -a \
        -o {output}
        """

use rule sort_and_convert as process_alignment with:
    output:
       "15_realignment/04_mapping/sorted/{sample}_R1R2_sorted_filtered_selected_accessions_sorted.bam",
    input:    
        rules.realign_reads.output,
    benchmark:
        "benchmarks/15_realignment/04_mapping/sorted/{sample}_R1R2_sorted_filtered_selected_accessions_sorted_benchmark.txt"

use rule index_bam as index_second_alignment with:
    output:
        "15_realignment/04_mapping/sorted/{sample}_R1R2_sorted_filtered_selected_accessions_sorted.bam.bai",
    input:    
        rules.process_alignment.output,
