from wordcloud import WordCloud

def generate_word_cloud(infile, outfile):
    """
      Return a word cloud based on taxon abundance data.
      Takes abundance data to generate an image illustrating taxon frequency
      in a sample.

      Parameters
      ----------
      infile : str
          Path where the abundance file file is located.
          File format either
            * 3 columns (e.g. species-filtered micop output):
                - first column: irrelevant (e.g. taxid)
                - second column: species name
                - third column: abundance/number of reads
            OR
            * 5 columns (e.g. consolidated and filtered species table):
                - first three columns same as above
                - fourth column: indication if secondary alignment was passed (irrelevant)
                - fifth column; indication if in addition coverage filter was passed
      outfile : str
          Path where the word clous should be saved in SVG format

      Returns
      -------
      N/A
      """
    # Create text from abundance dictionary holding species name as keys and
    # taxon abundances as values
    name2abundance = {}
    with open(infile, "r") as identifcation_table_in_fh:
        for line in identifcation_table_in_fh:
            columns = line.split(sep="\t")
            taxon_name, taxon_abundance = columns[5], columns[2]

            blast_identification = ""
            # If a 6-column input file is used
            if len(columns) > 3:
                blast_identification = columns[7].strip()

            # Only show species if a 6-columns input file is used and blast results
            # were passed or if a 3-column input file was used
            if blast_identification == "True" or len(columns) == 3:
                name2abundance[taxon_name.replace(" ", "-")] = int(taxon_abundance.strip())
            
        # If no species are identified set message to be displayed
        if len(name2abundance.keys()) == 0:
            name2abundance["no species identified"] = 1000

        # Generate word cloud based on word frequencies
        wordcloud = WordCloud(
            background_color="white",
            collocations=False,
            colormap="viridis",
            margin=0).fit_words(name2abundance)
        svg = wordcloud.to_svg()

        # Change display size in the browser
        svg = svg.replace('width="400" height="200"',
                    'width="100%" height="100%" viewBox="0 0 550 200"')
        with open(outfile, "w") as out_fh:
            out_fh.write(svg)

# identification_table_path = r"/data/projects/brama/pipeline_tests/validation/v0.26.2_v4.0.0_ATCC_v100-3M_005subVir_010subMicr_0002subHum/16_revalidation_blast/final_classification/simulated_reads_R1R2_sorted_filtered_selected_accessions_sorted_consensus_blast_final_identification.out"
# wordcloud_outpath = r"/data/projects/brama/pipeline_tests/validation/v0.26.2_v4.0.0_ATCC_v100-3M_005subVir_010subMicr_0002subHum/10_wordclouds/simulated_reads.html"
identification_table_path = snakemake.input[0]
wordcloud_outpath = snakemake.output[0]
generate_word_cloud(identification_table_path, wordcloud_outpath)