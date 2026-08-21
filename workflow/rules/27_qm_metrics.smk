# ToDo: change conda env to python 39
rule compute_quality_metrics:
    input:
        rules.create_final_identification_blast.output,
    output:
        report("quality_metrics/{sample}.jpg",
             category="quality control"),
        "quality_metrics/{sample}.txt",
    benchmark:
        "benchmarks/quality_metrics/{sample}_benchmark.txt"
    params:
        config["parameters"]["min_number_regions"],
        config["parameters"]["min_covbases"],
        config["parameters"]["min_region_length"],
        dynamic_path.databases.ground_truth,
    conda: "sm_python39",
    resources:
        mem_mb=300,
    script:
        "../scripts/qm_metrics.py"