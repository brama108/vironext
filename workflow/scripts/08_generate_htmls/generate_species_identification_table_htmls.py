#!/usr/bin/env python
# coding: utf-8

"""
Species Identification HTML Generator
=====================================

This script processes species identification data to generate a styled HTML table that summarizes
key taxonomic details and pipeline results. It takes as input a final identification table
and filtered BLAST results, filters the data based on reclassification and coverage criteria,
and creates an interactive HTML output.

The output HTML file includes:
- Detected taxa with links to NCBI taxonomy information.
- Visual indicators for reclassification (kaiju), coverage, and BLAST filtering outcomes.

Modules Used
------------
- `pandas`: For handling tabular data (optional, minimal use in this script).

Features
--------
1. Reads species identification data and BLAST filtering results.
2. Formats results into a visually appealing and interactive HTML table.
3. Highlights taxa passing or failing key filtering steps.
4. Includes links to NCBI taxonomy browser for quick reference.

Requirements
------------
- Input files:
  * A species identification table (TSV format).
  * Filtered BLAST results file (TSV format).

Parameters
----------
This script uses the following parameters via a workflow tool like Snakemake:
1. `species_identification_table_in_path` (str): Path to the input species identification table.
2. `species_identification_table_out_path` (str): Path to save the resulting HTML table.
"""

import pandas as pd

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
    
def create_species_table_html(
    species_table_file_in_path: str,
    species_coverage_table_file_in_path: str,
    unambiguous_reads_percentage_file_in_path: str,
    modified_species_table_file_out_path: str,
    html_head: str,
    html_tail: str,
    ):
    """
    Create a HTML from species-filtered MiCoP output.

    Parameters
    ----------
    species_table_file_in_path: str
        location of the final identfication table tsv file
    mpoified_species_table_file_out_path: str
        location of the final identification table HTML file
    species_coverage_table_file_in_path: str
        location of the coverage table needed to add %coverage to the HTML table
    html_head: str
        beginning of the HTML file
    html_tail: str
        end of the HTML file 
    """
    info_text = 'The following table lists the detected virus candidates. Taxa selected by the pipeline <span style="color: green; font-style: italic;">PASSED</span> reclassification as well as filtering.'
    species_coverage_table_df = pd.read_csv(species_coverage_table_file_in_path, sep="\t", names=["taxid", "coverage_percent"])
    unambiguous_reads_perc_df = pd.read_csv(unambiguous_reads_percentage_file_in_path, sep=" ", names=["taxid", "accession", "percentage_unambiguous_reads"])
    modified_species_table_file = open(modified_species_table_file_out_path, "a")
    modified_species_table_file.write(html_head)
    modified_species_table_file.write(f"""<h1>Sample: {snakemake.wildcards.sample}</h1><p>{info_text}</p>""")
    with open(species_table_file_in_path, "r") as species_table_fh:
        modified_species_table_file.write(f"""
    <table>
    <table style="border-collapse: collapse; width: 100%;">
    <tr style="background-color: #f2f2f2;">
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">TaxID</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Species name</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Blast hit</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Blast species</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Read number</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Reclassification</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Filtering</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Consolidation</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">Coverage</th>
    <th style="border: 1px solid #dddddd; text-align: left; padding: 8px;">High quality read assignment</th>
  </tr>
  """)
        for line in species_table_fh:
            values = line.split("\t")
            taxid = values[0]
            species_name = values[1]
            read_number = int(values[2])
            blast_classification = values[5]
            blast_hit = values[6] # this could be a taxon below species level
            # read_number = int(line.split()[-3])

            # extract percentage of unamiguously mapping reads for species taxid and assign color
            ambiguous_percent = unambiguous_reads_perc_df.query(f'taxid == {taxid}')["percentage_unambiguous_reads"].values[0] if taxid in unambiguous_reads_perc_df["taxid"].astype(str).values else "-"
            if ambiguous_percent == "-":
                ambiguous_percent_color = f'black">{ambiguous_percent}'
            else:
                if float(ambiguous_percent) < 15:
                    ambiguous_percent_color = f'red">{ambiguous_percent}%'
                elif float(ambiguous_percent) >= 15 and float(ambiguous_percent) < 20:
                    ambiguous_percent_color = f'orange">{ambiguous_percent}%'
                else:
                    ambiguous_percent_color = f'green">{ambiguous_percent}%'
            reclassification = 'green">PASSED' if values[3].strip()=="True" else 'red">FAILED'
            coverage_filter = 'green">PASSED' if values[4].strip()=="True" else 'red">FAILED'
            blast_filter = 'green">PASSED' if values[7].strip()=="True" else 'red">FAILED'
            #all_passed = True if "PASSED" in reclassification and "PASSED" in coverage_filter and "PASSED" in blast_filter else False
            
            # Extract coverage percent for taxid
            try:
                coverage_percent = f'{round(species_coverage_table_df[species_coverage_table_df["taxid"] == int(taxid)]["coverage_percent"].values[0],2)}%'
            except (IndexError, ValueError):
                coverage_percent = "-"  # if taxid not shown in the coverage plots (IndexError), if taxid=unclassified (ValueError)
            # Set link to NCBI
            link = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={taxid}&lvl=3&lin=f&keep=1&srchmode=1&unlock"

            # Set style for each row of the table
            row_style = "font-style: italic;" if values[7].strip()=="True" else ""

            # Write table to file and add thousand separator for read numbers (last columns)
            ## taxids with insufficient reads will not be schown in this table and all html files
            modified_species_table_file.write(f"""
    <tr style="{row_style}">
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px;"><a href={link}>{taxid}</a></td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px;">{species_name}</td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px;">{blast_hit}</td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px;">{blast_classification}</td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px;">{read_number:,}</td> 
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px; color: {reclassification}</td> 
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px; color: {coverage_filter}</td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px; color: {blast_filter}</td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px;">{coverage_percent}</td>
    <td style="border: 1px solid #dddddd; text-align: left; padding: 8px; color: {ambiguous_percent_color}</td>
  </tr>
  """)
            # else:
            #    exclude_taxids.append(taxid) 

        modified_species_table_file.write("</table>")
        modified_species_table_file.write(html_tail)
    modified_species_table_file.close()
    # with open(species_taxids_out_path, "w") as include_taxids_fout:
    #     include_taxids_fout.write("\n".join(include_taxids))

# basedir=r"/projects/brama/bloodvir/240308/v0.30.3"
# species_identification_table_in_path = rf"{basedir}/16_revalidation_blast/final_classification/AIL-138_2_S13_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification.out"
# species_coverage_table_in_path = rf"{basedir}/15_realignment/05_coverages/AIL-138_2_S13_coverage_percent.txt"
# percentage_unambiguous_reads = rf"{basedir}/15_realignment/06_ambiguous_alignment_stats/AIL-045_6_S11/AIL-045_6_S11_pct.txt"
# species_identification_table_out_path = rf"{basedir}/table.out.html"
species_identification_table_in_path = snakemake.input[0]
species_coverage_table_in_path = snakemake.input[1]
percentage_unambiguous_reads = snakemake.input[2]
species_identification_table_out_path = snakemake.output[0]

create_species_table_html(species_identification_table_in_path, species_coverage_table_in_path, percentage_unambiguous_reads, species_identification_table_out_path, HTML_HEAD, HTML_TAIL)



