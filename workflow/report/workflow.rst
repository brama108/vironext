Workflow information
====================
This workflow identifies species of a blood plasma virome that can be viewed in interactive reports. Adapters were removed with Cutadapt_.
Reads were assigned with `BWA MEM`_. Abundance was estimated using Microbial Community Profiler (MiCoP_). 
A graphical overview of the results can be found under `relative abundance`_, absolute values under `identified taxon list`_.
Information on coverages are stored in `coverages`_.

Quality control metrics are shown in `quality control`_.

.. _BWA MEM: http://bio-bwa.sourceforge.net
.. _Cutadapt: https://cutadapt.readthedocs.io
.. _MiCoP: https://github.com/smangul1/MiCoP


Filtering criteria
====================
Classifications are retained if at least one of the following filtering criteria is met: 

- Minimum number of non-adjacent regions: {{ snakemake.config["parameters"]["min_number_regions"] }}

OR

- Minimum number of covered bases: {{ snakemake.config["parameters"]["min_covbases"] }}

AND

- Minimum alignment length per region: {{ snakemake.config["parameters"]["min_region_length"] }}

Versioning
====================
- Pipeline version: {{ snakemake.config["pipeline_version"] }}
- Nucleotide database bwa: {{ snakemake.config["databases"]["RVDB_description"] }}
- Nucleotide database Blast: {{ snakemake.config["databases"]["RVDB_blast_description"] }}
- Protein database kaiju: {{ snakemake.config["protein_database_kaiju"] }}
- Protein database diamond: {{ snakemake.config["protein_database_diamond"] }}



