#!/usr/bin/env python
# Author: Martin Machyna
# Description: Script for filterg RVDB .fasta file based on mash distance

import argparse
import os
from subprocess import Popen, PIPE


parser = argparse.ArgumentParser(description='Script for filterg RVDB .fasta file based on mash distance')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-f', '--fasta', type=str, required=True, dest='fasta_file',
                    help='Annotated fasta file with headers in format: >rec|accesion_number|sequence_length|taxon_id|taxon_name')
requiredNamed.add_argument('-o', '--faOutDir', type=str, required=True, dest='fa_out_dir',
                    help='Directory for output .fasta files')
parser.add_argument('--filterLength', action='store_true', default = False, dest='filter_length',
                    help='Removes seqeunces < 1000 nt (<2000 nt for HIV-1)')





args = parser.parse_args()



def parse_fasta_taxons(fasta_file):
    """
    Parses annotated RVDB fasta file and groups sequences by taxon id
    Sequence records shorter than 1000 nt (shorter than 2000 nt for HIV)
    This will remove around 4,845 species with very short sequences and 409,159 HIV seqeunces
    input:
        fasta_file    : file annotate fasta file in format: >rec|accesion_number|sequence_length|taxon_id|taxon_name'
    output:
        taxon_dict    : dict {taxon_id: [[accession_number, sequence_length], ...]}
    """
    taxon_dict = {}

    fi = open(fasta_file, 'r')
    for line in fi:
        if line.startswith('>'):
            la = line.strip().split('|')

            acc = la[1]
            seq_len = la[2]
            tax_id = la[3]


            # Group sequnce records under taxon IDs
            if tax_id not in taxon_dict.keys():
                taxon_dict[tax_id] = [[acc, seq_len]]
            else:
                taxon_dict[tax_id].append([acc, seq_len])

    fi.close()

    return taxon_dict


def split_fasta_per_taxon(fasta_file, taxon_dict, filter_length = False):
    """
    Parses annotated RVDB fasta file and groups sequences by taxon id
    Sequence records shorter than 1000 nt (shorter than 2000 nt for HIV)
    This will remove around 4,845 species with very short sequences and 409,159 HIV seqeunces
    input:
        fasta_file    : file annotate fasta file in format: >rec|accesion_number|sequence_length|taxon_id|taxon_name'
        taxon_dict    : dict {taxon_id: [[accession_number, sequence_length], ...]}
        filter_length : boolean; Whether sequences shorter than 1000 nt (HIV 2000 nt) should be removed
    output:
        taxon_dict    : dict {taxon_id: [[accession_number, sequence_length], ...]}
    """
    taxon_info_dict = {}

    # Collect information about total number of records and number of RefSeq sequences
    for tax_id, records in taxon_dict.items():
        record_num = len(records)
        refseq_records = [rec[0] for rec in records if '_' in rec[0]]
        genome_fragments = len(refseq_records) # Species with fragmented genomes will have multiple refseq records

        taxon_info_dict[tax_id] = [record_num, genome_fragments]


    tempdir = 'tmp'
    os.mkdir(tempdir)

    fi = open(fasta_file, 'r')
    for line in fi:
        if line.startswith('>'):
            la = line.strip().split('|')

            acc = la[1]
            seq_len = la[2]
            tax_id = la[3]
            write = True


            # Set filter length for current species
            if '_' in acc:
                len_filter = 0    # Fiter for if the sequence is a RefSeq sequnces
            elif taxon_info_dict[tax_id][1] > 1:
                len_filter = 800  # Filter for multifragment species
            elif taxon_info_dict[tax_id][0] > 100000:
                len_filter = 2000 # Filter for HIV-1 with lot of records
            else:
                len_filter = 1000


            # Filter sequences
            if filter_length and int(seq_len) < len_filter:
                write = False
                taxon_dict[tax_id].remove([acc, seq_len])
                if len(taxon_dict[tax_id]) == 0:
                    del taxon_dict[tax_id]
                continue

            # Write sequence to the respective .fasta file
            fo_path = tempdir + '/' + tax_id + '.fa'
            fo = open(fo_path, 'a')

        if write:
            fo.write(line)

    fi.close()
    fo.close()

    return taxon_dict

# https://stackoverflow.com/questions/15343447/bash-style-process-substitution-with-pythons-popen
# https://stackoverflow.com/questions/163542/how-do-i-pass-a-string-into-subprocess-popen-using-the-stdin-argument
# https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc11  Bio.SeqIO.to_dict()



def get_mash_distance(taxon_id):
    """
    Calls mash dist to calculate mash distance and returns a list of mash distances
    input:
        taxon_id
     output:
        dict {ref_accessionr-qry_accession: mash_distnace}
    """
    hash_dict = {}
    tempfile = 'tmp/' + taxon_id + '.fa'

    p = Popen(["mash", "dist", "-p", "36", "-i", tempfile, tempfile], stdout=PIPE)

    output = p.communicate()[0]

    for line in output.splitlines():
        la = line.decode("utf-8").strip().split('\t')

        reference_id = la[0].split('|')[1]
        query_id = la[1].split('|')[1]
        mash_dist = la[2]

        hash_dict[reference_id + '-' + query_id] = mash_dist

    return hash_dict




def filter_grouped_records(taxons_dict, taxon_id):
    """
    Filter sequences inside species groups by their mash distance
    input:
        taxons_dict : dict {taxon_id: [[accession_number, sequence_length], ...]}
        taxon_id    : str, id of taxon to evaluate
     output:
        list [accession_number, ...]
    """
    def takeSecond(elem):
        return int(elem[1])

    final_list = []
    taxon_list = taxons_dict[taxon_id]
    hash_dict = get_mash_distance(taxon_id)

    # Sort descending by sequence length
    taxon_list.sort(key=takeSecond, reverse=True)

    # Set initial sequence for comparison to RefSeq records. RefSeq IDs alwas have '_' in them.
    final_list = [rec[0] for rec in taxon_list if '_' in rec[0]]
    # If there are no RefSeq records then use longest seqeunce for that taxon
    if len(final_list) == 0:
        longest_seq_acc = taxon_list[0][0]
        final_list = [longest_seq_acc]


    # Iteratively compare mash distances between the sequences
    for genome_record in taxon_list:
        accession = genome_record[0]

        # Check if distance to any of the already selected records is less than 0.15
        is_similar = [float(hash_dict[i + '-' + accession]) < 0.15 for i in final_list]

        # If all distances are more than 0.15 then add to list
        if sum(is_similar) == 0:
            final_list.append(accession)

    return final_list



def write_fasta_file(accession_list, faOutDir, fastaFile):
    """
    Filters .fasta files and outputs annotated fasta records
    input:
        accession_list  :   list;   list of accession numbers
        faOutDir        :   string; output folder
        fastaFile       :   string; path to annotated fasta file in format: >rec|accesion_number|sequence_length|taxon_id|taxon_name'
    """

    filename_out = faOutDir

    with open(filename_out, 'w') as fo:
        with open(fastaFile, 'r') as fi:

            for line in fi:
                if line.startswith('>'):
                    accesion = line.split('|')[1]

                if accesion in accession_list:
                    fo.write(line)







if __name__ == '__main__':
    taxon_dict = parse_fasta_taxons(args.fasta_file)
    taxon_dict_filt = split_fasta_per_taxon(args.fasta_file, taxon_dict, filter_length = args.filter_length)
        # This leaves about 16,407 species

    final_list = []
    for taxon_id in taxon_dict_filt.keys():
        filtered = filter_grouped_records(taxon_dict_filt, taxon_id)
        final_list.extend(filtered)

    # After filtering by length and mash score, we reduce RVDB to 42204 sequences distributed in 16407 species
    # HIV-1 reduced to  115 seqeunce records

    write_fasta_file(final_list, args.fa_out_dir, args.fasta_file)

