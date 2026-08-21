#!/usr/bin/env python
'''
Author: Martin Machyna
Email : martin.machyna@pei.de
Date: 2025-02-12
Description: Filters parts of non-viral sequences from RVDB .fasta file
'''


import argparse


def parse_arguments():
    """
    Parse CLI parameters.
    """
    parser = argparse.ArgumentParser(
        prog = 'Non-viral Filter',
        description='Script for filtering putative non-viral  sequences from RVDB .fasta files')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inFile', type=str, required=True, dest='in_file',
                        help='Path to RVDB .fasta input file')
    requiredNamed.add_argument('-o', '--outFile', type=str, required=True, dest='out_file',
                        help='Path to output .fasta file')
    requiredNamed.add_argument('-a', '--annotation', type=str, required=True, dest='annotation_file',
                        help='Annotation .tab file with non-viral sequences from RVDB website')

    args = parser.parse_args()

    return args



def parse_annotation(annotation_file):
    """
    Parses non-viral annotation .tab file
    :Parameters:
        annotation_file: str
            Path to .tab file

    :Returns:
        annot_dict: dict
            each accession has a list start and end coordinates
            {accession: [(start, end), (start, end), ...]}
            {str: list[tuple(int, int)]}
    """

    annot_dict = {}

    with open(annotation_file, 'r') as fi:
        next(fi)
        for line in fi:
            line = line.strip().split('\t')
            accession = line[0]
            start = int(line[2])
            end = int(line[3])

            if accession not in annot_dict:
                annot_dict[accession] = [(start,end)]
            else:
                annot_dict[accession].append((start,end))

    return annot_dict


def read_fasta(rvdb_fasta):
    """
    Parses multi-fasta file and stores it in dictionary
    :Parameters:
        rvdb_file: str
            Path to rvdb .fasta file
    :Returns:
        sequences: dict
            {header: sequence}
            {str: str}
    """
    sequences = {}

    with open(rvdb_fasta, 'r') as fi:
        for line in fi:
            if line.startswith('>'):
                header = line
                sequences[header] = ''
            else:
                sequences[header] += line.strip()

    return sequences


def mask(sequence, start, end, maskchar = '#'):
    """
    Masks part of string with masking character
    :Paramters:
        sequence: str
            String to mask
        start: int
            Start postion of interval to mask (0-based, closed)
        end: int
            End postition of interval to mask (0-based, open)
        maskchar: str
            Character used for masking

    :Returns:
        masked: str
            Masked string
    """
    if start not in range(len(sequence)) or end not in range(len(sequence)+1):
        raise ValueError('Masking coordinates are outside of string range')

    masked = sequence[:start] + (end-start)*maskchar + sequence[end:]

    return masked


def filter_fasta(sequences_dict, annot_dict, output_fasta):
    """
    Filters out segments of fasta sequences based on their start and end position
    :Parameters:
        sequences_dict : dict
            Dictionary storing fasta records
            {header: sequence}
            {str: str}
        annot_dict : dict
            Dictionary of non-viral segments annotation
            {accession: [(start, end), (start, end), ...]}
            {str: list[tuple(int, int)]}
        output_fasta : str
            Path for output .fasta file
    """

    fo = open(output_fasta, 'w')

    for header, sequence in sequences_dict.items():
        h_part = header.split('|')
        accession = h_part[2]

        if accession in annot_dict:
            # Mask all non-viral regions in given sequence
            for region in annot_dict[accession]:
                sequence = mask(sequence, region[0], region[1], '#')

            # Remove masked parts
            sequence = sequence.replace('#', '')

        if sequence != '':
            fo.write(header)
            fo.write(sequence + '\n')

    fo.close()



if __name__ == '__main__':
    args = parse_arguments()

    fasta_dict = read_fasta(args.in_file)
    annotation_dict = parse_annotation(args.annotation_file)

    filter_fasta(fasta_dict, annotation_dict, args.out_file)


