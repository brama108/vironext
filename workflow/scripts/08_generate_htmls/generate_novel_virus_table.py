#!/usr/bin/env python
# coding: utf-8

"""
@author: Martin Machyna
@date: Sep 26 2024
@desc: Script to convert .csv file into .html table
"""


import re
import sys
import os
import csv
import pandas as pd

def extract_taxon_name(x):
    "extract taxon name from subject title"
    if pd.notna(x):
        match = re.search(r"\[([^\]]+)\]", str(x))
        if match:
            return match.group(1)
    return x

HTML_HEAD = """
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
  </head>
  <body>
      <pre style="line-height: 1; font-family: monospace, monospace !important; font-size: 1em; white-space: pre; color: DarkBlue">
      """

HTML_TAIL = """
      </pre>
  </body>
</html>
"""

TABLE_HEAD = """
<table>
    <table style="border-collapse: collapse; width: 100%;">
    <tr style="background-color: #f2f2f2;">
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">id</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Contig name</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Length</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Viral frag</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Human frag</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Bact frag</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Viral / Total</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Classification (cross-check)</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Reads</th>
  </tr>
"""

# basedir="/projects/brama/pipeline_tests/validation/alilen150/v0.31.x/v0.31.1_v0.4.2_260213_5enteros_newBG"
# ML_results_path = f"{basedir}/13_ML_predict_virus/predicted_viruses_simulated_reads_BL_8.0/contigs_predicted.csv"
# blastx_result_path = f"{basedir}/13_ML_predict_virus/diamond/assembled_contigs_simulated_reads_BL_8.0/assembled_contigs_simulated_reads_BL_8.0_diamond.out"
# contig_remap_path = f"{basedir}/13_ML_predict_virus/contigs_remap/remap_stats_simulated_reads_BL_8.0.tsv"
# html_outpath = f"{basedir}/13_ML_predict_virus/report/simulated_reads_BL_8.0_novel_virus_table.html"
# read_cutoff = 100

ML_results_path = snakemake.input[0]
blastx_result_path = snakemake.input[1]
contig_remap_path = snakemake.input[2]
html_outpath = snakemake.output[0]
read_cutoff = snakemake.params[0]


# Check if all input files are present and non empty
for file_path in [ML_results_path, blastx_result_path, contig_remap_path]:
    if not (os.path.exists(file_path) and os.path.getsize(file_path) > 0):
        with open(html_outpath,'w') as html_outpath:
            html_outpath.write(HTML_HEAD)
            html_outpath.write('No novel viruses detected.<br>\nThere was likely insufficient number of unclassified reads to asseble genomic contigs.\n')
            html_outpath.write(HTML_TAIL)
        sys.exit(0)



# Load Ml novel virus prediction
ML_results_df = pd.read_csv(ML_results_path)

# load kraken2 results
blastx_df = pd.read_csv(blastx_result_path, header=None, usecols=[0, 2], names=['node_id', 'kingdom', 'subject_title'], sep='\t')

# Load read remapping results
contig_remap_df = pd.read_csv(contig_remap_path, sep='\t')

# Merge results
merged_df = pd.merge(ML_results_df, blastx_df, how="left", left_on="id", right_on="node_id")
merged_df = pd.merge(merged_df, contig_remap_df, how="left", left_on="id", right_on="#rname")
merged_df = merged_df[[
    'id',                 # contig id
    'length',
    '# viral fragments',
    '# plant fragments',
    '# bacterial fragments',
    '# viral / # total',  # viral / total fragments
    'decision',
    'numreads',             # Reads
    'meandepth',
    'subject_title',
]]

# filter out low depth contigs
merged_df = merged_df[(merged_df["numreads"] > read_cutoff) & (merged_df["decision"] == "virus")]

# extract taxon name from subject title
merged_df["taxon_name"] = merged_df["subject_title"].apply(extract_taxon_name)
merged_df.drop(columns=["meandepth", "subject_title", "decision"], inplace=True)
col = merged_df.pop("numreads")
merged_df.insert(7, "numreads", col)

# Open files for writing
html_outpath = open(html_outpath,'w')

# write header
html_outpath.write(HTML_HEAD)
html_outpath.write(f"""<h1>Sample: {snakemake.wildcards.sample}</h1>""")
html_outpath.write(TABLE_HEAD)

# generate table contents
for row in merged_df.itertuples():
    html_outpath.write(f'<tr style="border: 1px solid #dddddd; text-align: left; padding: 8px;">')

    # write table cells
    for col in row:
        html_outpath.write(f'<td>{col}</td>')
    
    html_outpath.write('</tr>')

# write end of the file
html_outpath.write('</table>')
html_outpath.write(HTML_TAIL)
html_outpath.close()


