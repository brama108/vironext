#!/usr/bin/env python
# Author: Martin Machyna
# Description: Script for processing RVDB .fasta file and annotating header with taxon id and length
#              Produces .fasta file with fasta header in fomrmat: >rec|accesion_number|sequence_length|taxon_id|taxon_name

import argparse


parser = argparse.ArgumentParser(description='Script for processing RVDB .fasta file and generating random viral metagenomic community')
parser.add_argument('fastaFiles', nargs='+', help='Fasta files to be labelled. Required.')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-a', '--acc2taxid', type=str, required=True, dest='acc2taxid_file',
                    help='Path to nucl_gb.accession2taxid file')
requiredNamed.add_argument('-t', '--names', type=str, required=True, dest='names_file',
                    help='Path to NCBI Taxonomy names.dmp file')
requiredNamed.add_argument('-n', '--nodes', type=str, required=True, dest='nodes_file',
                    help='Path to NCBI Taxonomy nodes.dmp file')
requiredNamed.add_argument('-o', '--faOutDir', type=str, required=True, dest='fa_out_dir',
                    help='Directory for output .fasta files')
parser.add_argument('-x', '--taxonLevel', type=str, default = 'species', dest='tax_level',
                    help='Taxon level at which taxons should be selected',
                    choices = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'])
parser.add_argument('--virusOnly', action='store_true', default = False, dest='virus_only',
                    help='Removes endogenous retroviruses and sequences not in viruses superkingdom')





args = parser.parse_args()


if not args.names_file or not args.acc2taxid_file or not args.fa_out_dir or not args.nodes_file:
    parser.print_help()
    exit(-1)

def parse_ncbi_taxonomy(acc2taxid_file, names_file, nodes_file):
    """
    Parses ncbi taxonomy dump files into dictionaries
    input:
        file nucl_gb.accession2taxid
        file names.dmp
        file nodes.dmp
     output:
        dict {accession: taxon_id}
        dict {taxon_id: scientific_name}
        dict {taxon_id: [parent_taxon_id, taxon_rank]}
    """

    acctaxid_dict, name_dict, node_dict = {}, {}, {}

    acc2taxid_fi = open(acc2taxid_file, 'r')
    for line in acc2taxid_fi:
        la = line.strip().split('\t')
        acc = la[0]
        taxid = la[2]
        acctaxid_dict[acc] = taxid
    acc2taxid_fi.close()

    names_fi = open(names_file, 'r')
    for line in names_fi:
        la = line.strip().split('\t')
        if la[6] == 'scientific name':
            node = la[0]
            sci_name = la[2]
            name_dict[node] = sci_name
    names_fi.close()

    nodes_fi = open(nodes_file, 'r')
    for line in nodes_fi:
        la = line.strip().split('\t')
        node = la[0]
        parent = la[2]
        rank = la[4]
        node_dict[node] = [parent, rank]
    nodes_fi.close()


    return acctaxid_dict, name_dict, node_dict



def find_path(node, node_dict, tax_level, virus_only):
    """
    Follows the taxonomic tree and makes a list of parent taxons
    Returns taxonID of species
    input:
        taxon id number to investigate
        dictionary {taxon_id: [paret_id, rank]}
        taxonomic level at which to report taxon id (one of 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
    output: taxon id
    """
    RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
    name_path = dict.fromkeys(RANKS, 'Unclassified')
    path = [node]
    while path[-1] != '1':
        if node_dict[path[-1]][1] in RANKS:
            name_path[node_dict[path[-1]][1]] = path[-1]
        path.append(node_dict[path[-1]][0])

    # Remove records that are not in virus taxonomic group
    if virus_only and name_path['superkingdom'] != '10239':
        name_path = dict.fromkeys(RANKS, 'not-virus')

    return name_path[tax_level]


def read_fasta(fasta_files):
    """
    Reads fasta file and produces dictionary of headres and sequence length
    input: list of paths to .fasta file from RVDB
    output: dictionary {record_number: [header, seqeunce_length]}
    """
    i = 0
    fasta_dict = {}

    for filename in fasta_files:
        with open(filename, 'r') as fi:

            for line in fi:
                line = line.strip()

                if line.startswith('>'):
                    i += 1
                    fasta_dict[i] = [line, 0]
                else:
                    fasta_dict[i][1] += len(line)

    return fasta_dict


def annotate_fasta(fasta_dict, acctaxid_dict, node_dict, name_dict, tax_level, virus_only):
    """
    Parses and annotates fasta records headers.
    input: dictionary {fasta_record_number: [header, seqeunce_length]}
           dictionary {accesion_number: taxon_id}
           dictionary {taxon_id: [parent_taxon_id, tax_tevel]}
           dictionary {taxon_id: scientific_name}
           string Taxonomy level at which the grouping should be performed ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
    output: dictionary {accesion_number: [fasta_record_number, accesion_number, sequence_length, taxon_id, taxon_name]}
    """
    annot_dict = {}

    for record_num in fasta_dict:
        header = fasta_dict[record_num][0].split('|') # >acc|NEIGHBOR|EU410304.1|Vaccinia virus GLV-1h68, complete genome.|Vaccinia virus GLV-1h68|VRL|29-SEP-2009
        seq_len = fasta_dict[record_num][1]

        seq_id = header[2]
        seq_id_base = seq_id.split('.')[0]

        try:
            record_tax_id = acctaxid_dict[seq_id_base]
            tax_id = find_path(record_tax_id, node_dict, tax_level, virus_only)
        except:
            print('Sequece ID ' + seq_id + ' not found in taxonomy dump file.')
            continue


        if tax_id == 'Unclassified':
            print('Filtered ', seq_id, ' record, because it did not have taxonID at ' + tax_level + ' level')
            continue
        elif tax_id == 'not-virus':
            print('Filtered ', seq_id, ' record, because it does not belong to Virus Superkingdom')
            continue


        tax_name = name_dict[tax_id]

        annot_dict[seq_id] = [record_num, seq_id, seq_len, tax_id, tax_name]

    return annot_dict


def write_fasta_file(annot_fasta_dict, faOutDir, fastaFiles):
    """
    Filters .fasta files and outputs annotated fasta records
    input:
        dictionary {accesion_number: [fasta_record_number, accesion_number, sequence_length, taxon_id, taxon_name]}
        output folder
        input rvdb.fasta file
    """

    filename_out = faOutDir

    with open(filename_out, 'w') as fo:
        for filename in fastaFiles:
            with open(filename, 'r') as fi:

                for line in fi:
                    if line.startswith('>'):
                        acc_num = line.split('|')[2] # >acc|NEIGHBOR|EU410304.1|Vaccinia virus GLV-1h68, complete genome.|Vaccinia virus GLV-1h68|VRL|29-SEP-2009


                        if acc_num in annot_fasta_dict.keys():
                            info = annot_fasta_dict[acc_num]
                            header = '>rec|' + info[1] + '|' + str(info[2]) + '|' + info[3] + '|' + info[4] + '\n'
                            fo.write(header)
                            out = True
                        else:
                            out = False

                    else:
                        if out:
                            fo.write(line)



if __name__ == '__main__':
    acctaxid_dict, name_dict, node_dict = parse_ncbi_taxonomy(args.acc2taxid_file, args.names_file, args.nodes_file)
    fasta_dict = read_fasta(args.fastaFiles)
    annot_fasta_dict = annotate_fasta(fasta_dict, acctaxid_dict, node_dict, name_dict, args.tax_level, args.virus_only)
    write_fasta_file(annot_fasta_dict, args.fa_out_dir, args.fastaFiles)
