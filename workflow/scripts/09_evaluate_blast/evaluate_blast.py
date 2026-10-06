#!/usr/bin/env python
# coding: utf-8

"""
Created on Apr 04 2023
@author: brama

BLAST Taxonomy Processor
========================

This script processes BLAST hit results, cross-references them with taxonomy 
data from an SQLite database, and generates an output file containing extended 
taxonomic information for each BLAST hit.

Purpose
-------
The script is designed to:
1. Load and process BLAST hit results.
2. Merge BLAST results with taxonomy data to enrich them with taxonomic information.
3. Identify and classify BLAST hits based on taxonomic congruence.
4. Write the processed and filtered results to an output file.

Requirements
------------
1. Input files:
   - SQLite database containing:
     * Taxonomic mapping: `taxonomy` table mapping TaxID to ScientificName.
     * Accession-to-TaxID mapping: `taxon_information_table`.
   - BLAST results file in tab-delimited format (output format 6).
2. Python packages:
   - `pandas`: For data manipulation.
   - `sqlite3`: For querying the SQLite database.

Input
-----
- SQLite Database Path:
  Path to the SQLite database containing taxonomy and accession-to-TaxID mappings.
- BLAST Hit Table:
  Path to the BLAST output file in tab-delimited format (output format 6).
  
Output
------
- Filtered and enriched BLAST hit results written to a tab-delimited file.

Classes
-------
1. BlastTaxonomyProcessor
   Encapsulates the functionality for processing BLAST hits and taxonomy data.
   - Methods:
     * `__init__`: Initializes the processor with database and file paths.
     * `load_taxonomy_data`: Loads taxonomy data into a dictionary and DataFrame.
     * `load_blast_hits`: Loads and filters the BLAST hit table.
     * `merge_taxonomy_and_blast`: Merges BLAST hits with taxonomy data and performs classification.
     * `write_identifications_to_file`: Writes the enriched and filtered data to a file.
"""

import sqlite3
import pandas as pd

class BlastTaxonomyProcessor:
    def __init__(self, sql_db_path, blast_hit_table_path, filtered_blast_out):
        """
        Initialize the processor with database and BLAST hit paths.
        
        Parameters:
        - sql_db_path (str): Path to the SQLite database.
        - blast_hit_table_path (str): Path to the BLAST hit table.
        - filtered_blast_out (str): Path to save the filtered BLAST output.
        """
        self.sql_db_path = sql_db_path
        self.blast_hit_table_path = blast_hit_table_path
        self.filtered_blast_out = filtered_blast_out  # Path for output
        self.len2taxid_merge_df = None
        self.blast_hits_df = None
        self.merged_df = None

    def load_taxonomy_data(self):
        """
        Load taxonomy data from the SQLite database.
        """
        conn = sqlite3.connect(self.sql_db_path)
        self.len2taxid_merge_df = pd.read_sql_query("SELECT * from taxon_information_table", conn)
        self.taxid2sciname_dict =  (
            pd.read_sql_query("SELECT * from taxonomy", conn)
            [["TaxID", "ScientificName"]].
            set_index("TaxID")["ScientificName"].
            to_dict()
        )
        conn.close()

    def load_blast_hits(self):
        """
        Load BLAST hits data and filter to retain the best hit per query.
        """
        self.blast_hits_df = pd.read_csv(
            self.blast_hit_table_path, sep="\t", names=[
                "query_name", "subject_name", "perc_identity", "alignment_length",
                "mismatches", "gap openings", "start of alignment in query",
                "end of alignment in query", "start of alignment in subject",
                "end of alignment in subject", "e-value", "bitscore", "query_len", "subject_len", "query_covs"
            ]
        )
        self.blast_hits_df = (
            self.blast_hits_df.sort_values(by=["bitscore", "e-value", "perc_identity"], 
                                           ascending=[False, True, False])
                               .groupby("query_name", as_index=False)
                               .first()
        )

    def merge_taxonomy_and_blast(self):
        """
        Merge BLAST hits with taxonomy data and add taxonomic information.
        """
        # Add taxid and species taxid of the query
        self.merged_df = (
            self.blast_hits_df.merge(
                self.len2taxid_merge_df[["Accession Version_DB", "taxid", "species_taxid"]],
                right_on="Accession Version_DB",
                left_on="query_name"
            )
            .drop(columns=["Accession Version_DB"])
            .rename(columns={"taxid": "q_taxid", "species_taxid": "q_species_taxid"})
        )

        # Add subject accession of the BLAST hit
        self.merged_df["subject_accession"] = self.merged_df["subject_name"].apply(
            lambda x: x.split("|")[2].split(".")[0]
        )

        # Add taxid and species taxid of the BLAST hit
        self.merged_df = (
            self.merged_df.merge(
                self.len2taxid_merge_df[["accession", "taxid", "species_taxid"]],
                right_on="accession",
                left_on="subject_accession"
            )
            .drop(columns=["accession"])
            .rename(columns={"taxid": "s_taxid", "species_taxid": "s_species_taxid"})
        )

        # add scientific name of blast hit species
        self.merged_df["blast_species_name"] = self.merged_df["s_species_taxid"].apply(lambda x: self.taxid2sciname_dict[x])
        
        # add scientific name of original blast hit 
        self.merged_df["blast_hit_name"] = self.merged_df["s_taxid"].apply(lambda x: self.taxid2sciname_dict[x])

    def write_identifications_to_file(self):
        """
        Write an extended blast hit list to a tsv file.
        """
        with open(self.filtered_blast_out, "w") as fout:
            self.merged_df.to_csv(fout, sep="\t", index=False)


# Main script logic
if __name__ == "__main__":
    # sql_db_path=r"/data/raw/Koenig/241223.bloodvir_MicopRequiredData/NCBI_taxonomy/80perc_contaRemoved_cond_C-RVDBv27_blast_dmpfiles.db" 
    # blast_hit_table_path=r"/data/projects/brama/pipeline_tests/validation/v0.25.0_v4.0.0_ATCC_v100-3M_010subVir_010subMicr_0002subHum/16_revalidation_blast/simulated_reads_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast.out"
    # filtered_blast_out_path = r"/data/projects/brama/pipeline_tests/validation/rvdbv29/test_kraken2/16_revalidation_blast/simulated_reads_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_filtered.out"
    blast_hit_table_path = snakemake.input[0]
    filtered_blast_out_path = snakemake.output[0]
    sql_db_path = snakemake.params[0]
    
    processor = BlastTaxonomyProcessor(sql_db_path, blast_hit_table_path, filtered_blast_out_path)
    processor.load_taxonomy_data()
    processor.load_blast_hits()
    processor.merge_taxonomy_and_blast()
    processor.write_identifications_to_file()
