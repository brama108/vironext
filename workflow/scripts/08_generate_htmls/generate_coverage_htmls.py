#!/usr/bin/env python
# coding: utf-8

"""
Created on Apr 04 2023

@author: brama
"""
"""
INTENDED USE
Create HTML file of the consolidated and coverage-filtered coverage plots.
Save mapping of taxids to coverage percentages in a text file.


DESCRIPTION
The text file is needed to add coverage information to the species information table.

REQUIREMENTS
This tool requires the following input files:
    * final samtools coverage plots

This script requires "pandas" to be installed within the Python environment the
script is executed in.

ASSOCIATED FUNCTIONS
N/A
"""

import os
import sqlite3
import re
import urllib.parse
import pandas as pd

def load_taxonomy(sql_db_path: str) -> tuple[pd.DataFrame, dict]:
    """
    Load taxonomy data from an SQLite database.

    This function reads two tables from the database:
    1. `taxonomy`: Contains the mapping of TaxIDs to their parent TaxIDs.
    2. `taxon_information_table`: Contains the mapping of Accessions to TaxIDs.

    Parameters
    ----------
    sql_db_path : str
        Path to the SQLite database file containing taxonomy information.

    Returns
    -------
    tuple[pd.DataFrame, dict]
        A tuple containing:
        - len2taxid_merge_df : pd.DataFrame
            DataFrame containing the contents of `taxon_information_table`.
        - node2parent_taxid_dict : dict
            Dictionary mapping TaxID to a list of parent TaxIDs.
    """
    conn = sqlite3.connect(sql_db_path)
    try:
        # Load taxonomy table into a DataFrame
        node2parent_taxid_df = pd.read_sql_query("SELECT * FROM taxonomy", conn)

        # Convert the taxonomy DataFrame into a dictionary with TaxID as the key
        node2parent_taxid_dict = (
            node2parent_taxid_df.set_index("TaxID")
            .agg(list, axis=1)
            .to_dict()
        )
    
        # Load taxon information table into a DataFrame
        len2taxid_merge_df = pd.read_sql_query("SELECT * FROM taxon_information_table", conn)
    finally:
        conn.close()

    return len2taxid_merge_df, node2parent_taxid_dict

def get_parent_taxid(taxid: str, node2parent_taxid_dict: dict) -> str:
                    """
                    Get the parent TaxID for a given TaxID.
                    """
                    return node2parent_taxid_dict.get(taxid, [None])[0]

HTML_HEAD = """
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        .header {
            line-height: 1;
            font-family: monospace, monospace !important;
            white-space: pre; 
            color: DarkBlue
        }
        .text {
            line-height: 1;
            font-family: monospace, monospace !important;
            font-size: 1em;
            white-space: pre;
            color: DarkBlue
        }
        .copy-btn {
            padding: 8px 16px;
            background-color: #007bff;
            color: white;
            border: none;
            cursor: pointer;
            font-size: 14px;
            margin-right: 10px;
            margin-bottom: 10px;
        }
        .accordion { }

        .copy-btn:hover {
            background-color: #0056b3;
        }

        .active, .accordion:hover {
            background-color: #ccc; 
            }

            .panel {
            padding: 0 18px;
            display: none;
            background-color: white;
            overflow: hidden;
        }
    </style>
  </head>
  <body>
      """

HTML_TAIL = """
    <script>
        function copyToClipboard(button) {
            var textToCopy = button.getAttribute("data-text");

            if (navigator.clipboard && navigator.clipboard.writeText) {
                // Use modern Clipboard API if available
                navigator.clipboard.writeText(textToCopy).catch(function(err) {
                    console.error("Failed to copy text: " + err);
                });
            } else {
                // Fallback for insecure context
                var textarea = document.createElement("textarea");
                textarea.value = textToCopy;
                document.body.appendChild(textarea);
                textarea.select();
                try {
                    document.execCommand("copy");
                } catch (err) {
                    console.error("Failed to copy text: " + err);
                }
                document.body.removeChild(textarea);
            }
        }

        function openBlast(url) {
            // Open the specified URL in a new tab
            window.open(url, '_blank');
        }
        var acc = document.getElementsByClassName("accordion");
        var i;

        for (i = 0; i < acc.length; i++) {
        acc[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var panel = this.nextElementSibling;
            if (panel.style.display === "block") {
            panel.style.display = "none";
            } else {
            panel.style.display = "block";
            }
        });
        }
    </script>
  </body>
</html>
"""

def is_empty(file_path: str) -> bool:
    """
    Check if a file is empty.

    Parameters
    ----------
    file_path : str
        The path to the file to check.

    Returns
    -------
    bool
        True if the file is empty, False otherwise.
    """
    if os.path.exists(file_path):
        if os.stat(file_path).st_size == 0:
            print("The file is empty.")
            return True
        else:
            return False
    else:
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")


def fasta2df(multifasta_path: str) -> pd.DataFrame:
    """
    Return a pd.DataFrame from a FASTA file.
    Converts a FASTA file to a DataFrame with header names and
        sequences.

    Parameters
    ----------
    multifasta_path : str
        location of the FASTA file containing consenss sequences
    Returns
    -------
    pd.DataFrame
        pd.DataFrame with sequence names as keys and DNA sequences
    """
    sequences = {}
    sequence_temp = ""
    with open(multifasta_path, "r", encoding="utf-8") as fin:
        for line in fin:
            if line.startswith(">"):
                try:
                    sequences[header] = sequence_temp
                    sequence_temp = ""
                    header = line.rstrip().strip(">")
                except NameError:
                    header = line.rstrip().strip(">")
            else:
                sequence_temp += line.rstrip().upper()
        sequences[header] = sequence_temp
    # Create dataframe from dictionary holding sequences and sequence headers
    db_df = pd.DataFrame(
        {
            "Header": sequences.keys(),
            "Sequence": sequences.values()
        }
        )
    return db_df

def create_coverage_graph_html(graph_file_in_path: str,
                                modified_graph_file_out_path: str,
                                coverage_percent_file_out_path: str,
                                html_head: str,
                                html_tail: str,
                                node2parent_taxid: dict,
                                len2taxid_merge: pd.DataFrame,
                                consensus_df: pd.DataFrame,
                                species_identification_table_in_path: str
                                ) -> None:
    """
    Create a HTML from a text file including species names.
    Samtools coverage --histogram / -m output only contains 
    accession number. This function adds the scientific name
    and creates an HTML file.
    Parameters
    ----------
    graph_file_in_path: str
        location of the samtools -m output
    modified_graph_file_out_path : str
        destination of the HTML file
    html_head: str
        beginning of the HTML file
    html_tail: str
        end of the HTML file 
    node2parent_taxid: dict
        dictionary holding the taxid as keys and parent taxids as values
    len2taxid_merge: pd.DataFrame
        Dataframe holding mapping between accessions and corresponding taxids
    Returns
    -------
    None
    """

    # load final species identification table
    final_identification_table_df = pd.read_csv(
        species_identification_table_in_path,
        sep="\t",
        names=[
            "taxid", "species", "read_count", "2nd_classification", "pass_cov_filter",
               "blast_identification", "blast_exact_hit", "pass_blast", "identity_perc",
               "alignment_length", "query_len", "subject_len", "query_covs", "accession"
               ],
        dtype={"taxid": str}
        )

    modified_graph_file = open(modified_graph_file_out_path, "a")
    modified_graph_file.write(html_head)
    modified_graph_file.write(f"<h1 class=header>Sample: {snakemake.wildcards.sample}</h1><p class=text>Coverage plots of the taxa found by the pipeline.</p><br><br>")
    with open(graph_file_in_path, "r") as graph_fh:
        i = 0
        for line in graph_fh:
            # Every 13th line of the samtools coverage output, there is an accession
            # followed by 11 lines of coverage information and 1 line with a space
            if i%13==0:
                # By default, taxa should not be skipped
                skip_lines = False

                # The accession number is separated by space from the sequence length
                accession = line.split()[0]
                sequence_length =line.split()[1]
                
                # Get TaxID
                to_be_appended_taxids = len2taxid_merge.query(f"`Accession Version_DB` == '{accession}'")["taxid"].values[0]
                
                # Get scientific name
                to_be_appended_name = node2parent_taxid[str(to_be_appended_taxids)][2]

                # Get parent TaxID (make sure to catch species level taxid)
                to_be_appended_taxids = [to_be_appended_taxids]
                to_be_appended_taxids.append(get_parent_taxid(str(to_be_appended_taxids[-1]), node2parent_taxid))

                # Get parent of parent TaxID (make sure to catch species level taxid)
                to_be_appended_taxids.append(get_parent_taxid(str(to_be_appended_taxids[-1]), node2parent_taxid))

                # Get parent of parent of parent TaxID (make sure to catch species level taxid)
                to_be_appended_taxids.append(get_parent_taxid(str(to_be_appended_taxids[-1]), node2parent_taxid))

                # Get the consensus sequence
                try:
                    consensus_sequence = consensus_df.query(f"Header == '{accession}'")['Sequence'].iloc[0]
                except IndexError:
                    consensus_sequence = "No consensus sequence was generated for this accession. This is likely due to insufficient data."

                # # Exclude taxa with insufficient reads from reporting
                # if not any(taxon in include_taxa for taxon in [to_be_appended_taxid, to_be_appended_parent_taxid, to_be_appended_parent_parent_taxid]):
                #     # Skip taxon with insufficient reads
                #     skip_lines = True
                #     i += 1
                #     continue
                
                # Set links to NCBI
                link_accession = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
                link_taxid = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={to_be_appended_taxids[0]}&lvl=3&lin=f&keep=1&srchmode=1&unlock"
                link_parent_taxid = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={to_be_appended_taxids[1]}&lvl=3&lin=f&keep=1&srchmode=1&unlock"
                link_parent_parent_taxid = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={to_be_appended_taxids[2]}&lvl=3&lin=f&keep=1&srchmode=1&unlock"
                link_parent_parent_parent_taxid = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={to_be_appended_taxids[3]}&lvl=3&lin=f&keep=1&srchmode=1&unlock"
                
                # Label species that passed the blast filter with "IDENTIFIED" and those that did not with "FILTERED OUT"
                final_identification_table_filtered_df = final_identification_table_df[final_identification_table_df["pass_blast"]==True]
                taxid_set = set(final_identification_table_filtered_df["taxid"].values)
                filters_passed = "IDENTIFIED" if any(t in taxid_set for t in to_be_appended_taxids) else "FILTERED OUT"
                # Form the strings to be appended to html file
                to_be_appended = f"<b>{filters_passed}</b> - {to_be_appended_name} (<a href={link_accession}>{accession}</a>, TaxID <a href={link_taxid}>{to_be_appended_taxids[0]}</a>, 1-level-up TaxID <a href={link_parent_taxid}>{to_be_appended_taxids[1]}</a>, 2-levels-up TaxID <a href={link_parent_parent_taxid}>{to_be_appended_taxids[2]}</a>, 3-levels-up TaxID <a href={link_parent_parent_parent_taxid}>{to_be_appended_taxids[3]}</a>) {sequence_length}\n"
                
                fasta_header = f">{snakemake.wildcards.sample} {accession} {to_be_appended_name} taxid {to_be_appended_taxids[0]}"
                fasta_sequence = f"{fasta_header}\n{consensus_sequence}"
                sequence_fasta_format = urllib.parse.quote(fasta_sequence)

                
                modified_graph_file.write("<pre class=text>")
                modified_graph_file.write(to_be_appended)
            else:
                if skip_lines == True:
                    i += 1
                    continue
                else:
                    to_be_appended = line
                    modified_graph_file.write(to_be_appended)
                    if (i + 2)%13==0:
                        modified_graph_file.write("</pre>")
                        # add button to copy consensus sequence
                        modified_graph_file.write(f'<button class="copy-btn" data-text="{fasta_sequence}" onclick="copyToClipboard(this)">Copy consensus sequence</button>')
                        modified_graph_file.write(
    f'<button class="copy-btn" onclick="openBlast(\'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome\')">Open Blast</button>'
)
                        modified_graph_file.write(
    f'<button class="copy-btn" onclick="openBlast(\'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&LAYOUT=OneWindows&AUTO_FORMAT=Fullauto&PAGE=Nucleotides&NCBI_GI=yes&FILTER=L&HITLIST_SIZE=100&SHOW_OVERVIEW=yes&AUTO_FORMAT=yes&SHOW_LINKOUT=yes&QUERY={sequence_fasta_format}\')">Blast consensus sequence</button>'
)
                        modified_graph_file.write("<p><u class=text>BLAST results:</u></p>")
                        blast_info = final_identification_table_df[final_identification_table_df['taxid'].isin(to_be_appended_taxids)]
                        link_blast_accession = f"https://www.ncbi.nlm.nih.gov/nuccore/{blast_info['accession'].values[0]}"
                        modified_graph_file.write(f"<span class=text><a href={link_blast_accession}>{blast_info['accession'].values[0]}</a> - {blast_info['blast_exact_hit'].values[0]}, Identity: {blast_info['identity_perc'].values[0]}%, Consensus sequence length: {blast_info['query_len'].values[0]} bp, Reference sequence length: {blast_info['subject_len'].values[0]} bp, Consensus sequence coverage: {blast_info['query_covs'].values[0]}% \n")              
                    # write coverage% and taxids to file as input for generation of specied identification table
                    elif (i+9)%13==0:
                         match = re.search(r'Percent covered:\s*([\d.]+)%', to_be_appended)
                         percent_covered = match.group(1)
                         with open(coverage_percent_file_out_path, "a") as fout:
                              fout.write("\n".join(f"{taxid}\t{percent_covered}" for taxid in to_be_appended_taxids) + "\n")

                    else:
                        pass
            i += 1
        modified_graph_file.write(html_tail)
    modified_graph_file.close()
# coverages_graph_file_path = r"/projects/brama/bloodvir/230100_NYC/batch4_alilen0/15_realignment/05_coverages/437_S12_L001_R1R2_sorted_filtered_selected_accessions_sorted.graph.txt"
# consensus_sequences_path = r"/projects/brama/bloodvir/230100_NYC/batch4_alilen0/15_realignment/06_consenus_sequence/437_S12_L001_R1R2_sorted_filtered_selected_accessions_sorted_consensus.fasta"
# final_species_identification_table_path = r"/projects/brama/bloodvir/230100_NYC/batch4_alilen0/16_revalidation_blast/final_classification/437_S12_L001_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification.out"
# coverag_plot_html_out_path = r"/projects/brama/bloodvir/230100_NYC/batch4_alilen0/15_realignment/05_coverages/437_S12_L001_2.html"
# coverage_percent_out_path = r"/projects/brama/bloodvir/230100_NYC/batch4_alilen0/15_realignment/05_coverages/coverage_percent.txt"
# taxonomy_database = r"/data/raw/Koenig/250128.bloodvir_MicopRequiredData/NCBI_taxonomy/contaRemoved_condC-RVDBv29.0_dmpfiles.db"
coverages_graph_file_path = snakemake.input[0]
consensus_sequences_path = snakemake.input[1]
final_species_identification_table_path = snakemake.input[2]
coverag_plot_html_out_path = snakemake.output[0]
coverage_percent_out_path = snakemake.output[1]
taxonomy_database = snakemake.params[0]

if not is_empty(coverages_graph_file_path):
    consensus_seqs_df = fasta2df(consensus_sequences_path)
    len2taxid_df, node2par_taxid_dict = load_taxonomy(taxonomy_database)
    create_coverage_graph_html(coverages_graph_file_path, coverag_plot_html_out_path, coverage_percent_out_path, HTML_HEAD, HTML_TAIL, node2par_taxid_dict, len2taxid_df, consensus_seqs_df, final_species_identification_table_path)
else:
    with open(coverag_plot_html_out_path, "w") as coverage_plots_html_fout:
        coverage_plots_html_fout.write("no plots to show")
    with open(coverage_percent_out_path, "w") as coverage_percentages_fout:
        coverage_percentages_fout.write("no coverages to be computed")
