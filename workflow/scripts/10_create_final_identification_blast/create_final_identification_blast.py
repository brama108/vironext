#!/usr/bin/env python
# coding: utf-8

"""
Created on Apr 04 2023

@author: brama
"""
"""
INTENDED USE
Filter identification table after bwa/micop classification using filtered blast results.

...
"""
import pandas as pd

def add_blast_information(blast_information_type: str) -> None:
    """
    Add blast information to the final identification table.
    Parameters
    ----------
    blast_information_type: str
        May be one of the following:
        - "blast_species_name"
        - "blast_hit_name"
        - "perc_identity"
        - "alignment_length"
        - "query_len"
        - "subject_len"
        - "query_covs"
        - "accession"

    Return
    None
    ----------
    """
    if blast_information_type == "blast_perc_identity":
        blast_processed_results_df[blast_information_type] = round(blast_processed_results_df[blast_information_type].values[0], 1)
    # add blast information to final identification table
    ## if statement necessary to avoid crash when blast search has not been done (i.e. taxid is empty)
    final_identification_table_df[blast_information_type] = final_identification_table_df["taxid"].apply(
    lambda x: blast_processed_results_df[blast_processed_results_df["q_species_taxid"] == x][blast_information_type].values[0]
    if not blast_processed_results_df[blast_processed_results_df["q_species_taxid"] == x].empty
    else "-"
)

path_to_filtered_micop_kaiju_results = snakemake.input[0]
blast_filtered_results_path = snakemake.input[1]
outpath_to_final_identification_table = snakemake.output[0]
# path_to_filtered_micop_kaiju_results = r"/data/projects/brama/pipeline_tests/validation/v0.25.0_v4.0.0_ATCC_v100-3M_010subVir_010subMicr_0002subHum/14_validate_calls/filtered_identifications/simulated_reads_R1R2_sorted_filtered_abundances_regions_micop_kaiju_filtered.txt"
# blast_filtered_results_path = r"/data/projects/brama/pipeline_tests/validation/v0.25.0_v4.0.0_ATCC_v100-3M_010subVir_010subMicr_0002subHum/16_revalidation_blast/simulated_reads_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_filtered.out"
# outpath_to_final_identification_table = r"/data/projects/brama/pipeline_tests/validation/v0.25.0_v4.0.0_ATCC_v100-3M_010subVir_010subMicr_0002subHum/16_revalidation_blast/final_classification/simulated_reads_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification2.out"


# load filtered blast output
blast_processed_results_df = pd.read_csv(blast_filtered_results_path, sep="\t", dtype={"q_species_taxid": str})

# add blast accession columns to blast processed results
blast_processed_results_df["accession"] = blast_processed_results_df["subject_name"].str.split("|").str[2]

# load consolidated micop/kaiju results
filtered_micop_kaiju_results_df = pd.read_csv(path_to_filtered_micop_kaiju_results, sep = "\t", dtype={"taxid": str, "species": str, "read_count": int}, names = ["taxid", "species", "read_count", "2nd_classification", "pass_cov_filter"])

final_identification_table_df = filtered_micop_kaiju_results_df.copy(deep=True)

# extract taxids of blast classification
include_taxids_blast = list(map(str, blast_processed_results_df["s_species_taxid"]))

add_blast_information("blast_species_name")
add_blast_information("blast_hit_name")

# Indicate if species idendified by blast is the same as identified by micop/kaiju
final_identification_table_df["blast"] = final_identification_table_df.apply(
    lambda x: True if x["species"] == x["blast_species_name"] else False, axis=1
)

add_blast_information("perc_identity")
add_blast_information("alignment_length")
add_blast_information("query_len")
add_blast_information("subject_len")
add_blast_information("query_covs")
add_blast_information("accession")


final_identification_table_df.sort_values(["blast","pass_cov_filter", "2nd_classification", "read_count"], ascending=False, inplace=True)

# write results to file
final_identification_table_df.to_csv(outpath_to_final_identification_table, header=False, sep="\t", index=False)
