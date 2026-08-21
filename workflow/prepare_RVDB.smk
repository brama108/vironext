configfile: workflow.source_path("../config/config_prepare_RVDB.yaml")

rule all:
    input:
        "db_out/C-RVDBv" + config["RVDB"] + "_annotated_filtered.fasta"

rule get_taxonomy:
    output:
        nodes = "db_in/nodes.dmp",
        names = "db_in/names.dmp"
    threads: 1
    shell:
        """
        wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O db_in/taxdump.tar.gz
        wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5 -O db_in/taxdump.tar.gz.md5

        echo "$(cut -f1 -d ' ' db_in/taxdump.tar.gz.md5)  db_in/taxdump.tar.gz" | md5sum -c -

        tar -xvzf db_in/taxdump.tar.gz -C db_in/
        """

rule get_accessions:
    output:
        acc2tax = "db_in/nucl_gb.accession2taxid"
    threads: 1
    shell:
        """
        wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz -O db_in/nucl_gb.accession2taxid.gz
        wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5 -O db_in/nucl_gb.accession2taxid.gz.md5

        echo "$(cut -f1 -d ' ' db_in/nucl_gb.accession2taxid.gz.md5)  db_in/nucl_gb.accession2taxid.gz" | md5sum -c -

        gzip -d db_in/nucl_gb.accession2taxid.gz
        """

rule get_rvdb:
    output:
        output_file = "db_in/C-RVDBv" + config["RVDB"] + ".fasta"
    params:
        rvdb_version = config["RVDB"]
    threads: 1
    shell:
        """
        wget -q https://rvdb.dbi.udel.edu/download/C-RVDBv{params.rvdb_version}.fasta.gz -O db_in/C-RVDBv{params.rvdb_version}.fasta.gz

            #TODO: add md5 sum check step

        gzip -d {output.output_file}.gz
        """
rule get_rvdb_nonviral:
    output:
        nonviral_annotation = "db_in/putative_non_viral_annotation.tab"
    threads: 1
    shell:
        """
        file_url=$(wget -q -O - https://rvdb.dbi.udel.edu/previous-release | grep 'far fa-file-alt' -B 1 | grep '<a' | head -n 1 | cut -d '"' -f 4)
        #Note: This is a workaround until they make file names consistent. It might brake if they change their website code.

        wget -q -O db_in/putative_non_viral_annotation.tab ${{file_url}}
        """

rule filter_nonviral:
    input:
        rvdb = lambda wildcards: config["RVDB_fasta"] if config["RVDB_fasta"] != "" else rules.get_rvdb.output.output_file,
        nv_annot = rules.get_rvdb_nonviral.output.nonviral_annotation
    output:
        rvdb_nv = "db_out/C-RVDBv" + config["RVDB"] + "_fnv.fasta"
    params:
        filter_nv = workflow.source_path("scripts/parse_rvdb/filter_nonviral.py")
    threads: 1
    shell:
        """
        python {params.filter_nv} \
            -i {input.rvdb} \
            -o {output.rvdb_nv} \
            -a {input.nv_annot}
        """


rule annotate_rvdb:
    input:
        rvdb = rules.filter_nonviral.output.rvdb_nv,
        acc2taxid = lambda wildcards: config["acc2taxid"] if config["acc2taxid"] != "" else rules.get_accessions.output.acc2tax,
        taxNodes = lambda wildcards: config["nodes"] if config["nodes"] != "" else rules.get_taxonomy.output.nodes,
        taxNames = lambda wildcards: config["names"] if config["names"] != "" else rules.get_taxonomy.output.names
    output:
        "db_out/C-RVDBv" + config["RVDB"] + "_annotated.fasta"
    params:
        annotate_rvdb=workflow.source_path("scripts/parse_rvdb/annotate_rvdb.py")
    threads: 1
    shell:
        """
        python {params.annotate_rvdb} {input.rvdb} \
            -a {input.acc2taxid} \
            -t {input.taxNames} \
            -n {input.taxNodes} \
            -o {output} \
            --taxonLevel species \
            --virusOnly
        """

rule filter_rvdb:
    input:
        "db_out/C-RVDBv" + config["RVDB"] + "_annotated.fasta"
    output:
        "db_out/C-RVDBv" + config["RVDB"] + "_annotated_filtered.fasta"
    params:
        filter_rvdb=workflow.source_path("scripts/parse_rvdb/filter_rvdb.py")
    threads: 30
    resources:
        mem_mb=900000
    conda:
        "sm_mash"
    shell:
        """
        python {params.filter_rvdb} -f {input} \
                     -o {output} \
                     --filterLength

        rm -r tmp/
        """



