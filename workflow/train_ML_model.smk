configfile: workflow.source_path("../config/config_train_ML_model.yaml")


rule all:
    input:
        "model_weights/model_5_500.h5",
        "model_weights/model_5_1000.h5",
        "model_weights/model_7_500.h5",
        "model_weights/model_7_1000.h5",
        "model_weights/model_10_500.h5",
        "model_weights/model_10_1000.h5"

rule get_human_genome:
    output:
        output_file = "db_in/chm13v2.0.fasta"
    threads: 1
    shell:
        """
        wget -q https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -O db_in/chm13v2.0.fa.gz

        echo "06caded0055ec647a69acb708e13beff db_in/chm13v2.0.fa.gz" | md5sum -c -

        gzip -d db_in/chm13v2.0.fa.gz
        mv db_in/chm13v2.0.fa db_in/chm13v2.0.fasta
        """

rule get_bacteria_genomes:
    output:
        output_file = "db_in/bacteria.fasta"
    threads: 30
    conda:
        "sm_parallel"
    shell:
        """
        wget -O- https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt \
            | awk -v "FS=\t" '$11 == "latest" && $12 == "Complete Genome" {{print $20}}' \
            > db_in/bact_url_download.txt

        parallel -j 36 wget -qO- {{}}/{{/}}_genomic.fna.gz :::: db_in/bact_url_download.txt > db_in/bacteria.fasta.gz

        gzip -d db_in/bacteria.fasta.gz
        """

rule ML_prepare_data:
    input:
        virus_fa = config["RVDB_fasta"],
        human_fa = lambda wildcards: config["human_fasta"] if config["human_fasta"] != "" else rules.get_human_genome.output.output_file,
        bacteria_fa = lambda wildcards: config["bacteria_fasta"] if config["bacteria_fasta"] != "" else rules.get_bacteria_genomes.output.output_file
    output:
        "train_data/encoded_train_500.hdf5",
        "train_data/encoded_train_1000.hdf5"
    params:
        vh_prepare=workflow.source_path("scripts/virhunter/prepare_ds.py"),
        vh_train=workflow.source_path("scripts/virhunter/train.py"),
        vh_utils=workflow.source_path("scripts/virhunter/utils/preprocess.py"),
        vh_mod10=workflow.source_path("scripts/virhunter/models/model_10.py"),
        vh_mod7=workflow.source_path("scripts/virhunter/models/model_7.py"),
        vh_mod5=workflow.source_path("scripts/virhunter/models/model_5.py"),
        vh_config=workflow.source_path("../config/config_virhunter.yml")
    conda:
        "sm_virhunter"
    threads: 60
    resources:
        mem_mb=50000
    shell: 
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/

        sed -e 's|<A>|{input.virus_fa}|' \
            -e 's|<B>|{input.human_fa}|' \
            -e 's|<C>|{input.bacteria_fa}|' \
            -e 's|<D>|train_data|' \
            -e 's|<E>|train_data|' \
            -e 's|<F>|model_weights|' \
            {params.vh_config} > train_data/config_virhunter.yml

        python {params.vh_prepare} train_data/config_virhunter.yml
        """

rule ML_train:
    input:
        "train_data/encoded_train_500.hdf5",
        "train_data/encoded_train_1000.hdf5"
    output:
        "model_weights/model_5_500.h5",
        "model_weights/model_5_1000.h5",
        "model_weights/model_7_500.h5",
        "model_weights/model_7_1000.h5",
        "model_weights/model_10_500.h5",
        "model_weights/model_10_1000.h5"
    params:
        vh_prepare=workflow.source_path("scripts/virhunter/prepare_ds.py"),
        vh_predict=workflow.source_path("scripts/virhunter/predict.py"),
        vh_train=workflow.source_path("scripts/virhunter/train.py"),
        vh_utils=workflow.source_path("scripts/virhunter/utils/preprocess.py"),
        vh_utils2=workflow.source_path("scripts/virhunter/utils/batch_loader.py"),
        vh_mod10=workflow.source_path("scripts/virhunter/models/model_10.py"),
        vh_mod7=workflow.source_path("scripts/virhunter/models/model_7.py"),
        vh_mod5=workflow.source_path("scripts/virhunter/models/model_5.py"),
        vh_config=workflow.source_path("../config/config_virhunter.yml")
    conda:
        "sm_virhunter"
    threads: 10
    resources:
        mem_mb=100000,
        partition="GPU"
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/

        python {params.vh_train} train_data/config_virhunter.yml
        """




