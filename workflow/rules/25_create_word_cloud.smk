rule create_word_cloud:
    """
    Generate wordcloud depicting the relative abundance of identified viruses
    """
    output:
        report("10_wordclouds/{sample}.html", caption="../report/word_cloud_caption.txt", category="relative abundance", subcategory="word cloud" ),
    input:    
        rules.create_final_identification_blast.output,
    benchmark:
        "benchmarks/10_wordclouds/{sample}_benchmark.txt"
    threads: 1
    resources:
        mem_mb=100
    conda: "sm_wordcloud"
    script:
        "../scripts/07_generate_wordcloud/generate_wordcloud.py"