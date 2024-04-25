#!/usr/bin/env python3
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process VCF file and save output')
    parser.add_argument('--tsv', type=str, help='Path to input tsv file', required=True)
    parser.add_argument('--infasta', type=str, help='Path to input fasta file', required=True)
    parser.add_argument('--outfasta', type=str, help='Path to modified fasta file', required=True)
    return parser.parse_args()


def replace_nucleotides_with_iupac(fasta_file, tsv_file):
    # Read the TSV file and extract the Position and IUPAC information
    variants = []
    with open(tsv_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            Chromosome,Position,Reference,Combined_Alternate,Combined_Alleles,IUPAC = line.strip().split('\t')
            variants.append((int(Position), IUPAC))

    # Read the FASTA file
    with open(fasta_file, 'r') as file:
        lines = file.readlines()

    # Iterate over the variants
    for variant in variants:
        Position = variant[0]
        IUPAC = variant[1]

        # Find the sequence line in the FASTA file
        sequence_line = None
        for line in lines:
            if not line.startswith('>'):  # Skip header lines
                sequence_line = lines.index(line)
                break

        if sequence_line is not None:
            sequence = list(lines[sequence_line].strip())

            # Replace the nucleotide at the specified Position with the IUPAC code
            sequence[Position - 1] = IUPAC

            # Update the sequence line in the FASTA file
            lines[sequence_line] = ''.join(sequence) + '\n'

    # Write the modified FASTA file
    modified_file = args.outfasta
    with open(modified_file, 'w') as file:
        file.writelines(lines)

    print("Modified FASTA file created:", modified_file)



if __name__ == '__main__':
    args = parse_arguments()
    tsv_file = args.tsv
    fasta_file = args.infasta
    replace_nucleotides_with_iupac(fasta_file, tsv_file)
