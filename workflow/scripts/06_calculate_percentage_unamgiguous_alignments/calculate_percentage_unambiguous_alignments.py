#!/usr/bin/env python
# coding: utf-8

"""
Created on Apr 04 2023

@author: brama
"""
"""
INTENDED USE
Calculate percentage of unambiguous alignments from total alignments.


DESCRIPTION
This script takes the total number of alignments and the number of unambiguous
alignments and calculates the percentage of unambiguous alignments.

REQUIREMENTS
This tool requires the following input files:
    * total number of alignments per spieces
    * number of unambiguous alignments per species

ASSOCIATED FUNCTIONS
This file contains the following functions:
    * N/A

"""
import os

# Snakemake inputs
path_to_total_number_of_alignments = snakemake.input.total
path_to_number_of_unambig_alignments = snakemake.input.unambig
path_to_percentage_output = snakemake.output[0]
accession_dir = snakemake.params.accession_dir  # path to directory listing the accession that have been aligned to 

# Open output
with open(path_to_percentage_output, 'w') as out:
    # Loop over all accessions
    for t_file, u_file in zip(path_to_total_number_of_alignments, path_to_number_of_unambig_alignments):
        accession = os.path.basename(t_file).replace("_total_number_alignments.txt","")
        
        # Map accession to species taxid
        accession_file = os.path.join(accession_dir, f"{accession}.txt")
        with open(accession_file) as f:
            species_taxid, acc_from_file = f.read().strip().split()
        
        # Read counts
        with open(t_file) as f:
            total_count = int(f.read().strip())
        with open(u_file) as f:
            unambig_count = int(f.read().strip())
        
        # Calculate percentage safely
        pct = (unambig_count / total_count * 100) if total_count > 0 else 0.0
        
        # Write species_taxid, accession, percentage
        out.write(f"{species_taxid} {accession} {pct:.2f}\n")
