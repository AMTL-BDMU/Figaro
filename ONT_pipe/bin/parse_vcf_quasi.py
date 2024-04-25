#!/usr/bin/env python3
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process VCF file and save output')
    parser.add_argument('--vcf', type=str, help='Path to input VCF file', required=True)
    parser.add_argument('--out', type=str, help='Path to output file', required=True)
    return parser.parse_args()


iupac_codes = {
    ('A', 'G'): 'R',
    ('C', 'T'): 'Y',
    ('G', 'C'): 'S',
    ('A', 'T'): 'W',
    ('G', 'T'): 'K',
    ('A', 'C'): 'M',
    ('C', 'G', 'T'): 'B',
    ('A', 'G', 'T'): 'D',
    ('A', 'C', 'T'): 'H',
    ('A', 'C', 'G'): 'V',
    ('A', 'C', 'G', 'T'): 'N',
    ('N', 'A', 'G'): 'R',
    ('N', 'C', 'T'): 'Y',
    ('N', 'G', 'C'): 'S',
    ('N', 'A', 'T'): 'W',
    ('N', 'G', 'T'): 'K',
    ('N', 'A', 'C'): 'M',
    ('N', 'C', 'G', 'T'): 'B',
    ('N', 'A', 'G', 'T'): 'D',
    ('N', 'A', 'C', 'T'): 'H',
    ('N', 'A', 'C', 'G'): 'V',
    ('A'): 'A',
    ('G'): 'G',
    ('C'): 'C',
    ('T'): 'T',
    ('N', 'A'): 'A',
    ('N', 'G'): 'G',
    ('N', 'C'): 'C',
    ('N', 'T'): 'T'
}

def find_iupac_code(combined_alleles):
    combined_alleles_set = set(combined_alleles)
    for alleles, iupac_code in iupac_codes.items():
        if set(alleles) == combined_alleles_set:
            return iupac_code
    return 'N'

def parse_vcf(file_path, output_file):
    variants = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                # Skip header lines
                continue
            fields = line.split('\t')
            filters = fields[7].split(';')
            af_field = next((field for field in filters if field.startswith('AF=')), None)
            
            # Check if "PASS" is present in the filters
            if af_field is not None:
                af = float(af_field.split('=')[1])
                if af >= 0.2:
                    chromosome = fields[0]
                    position = int(fields[1])
                    ref = fields[3].upper()
                    alt = fields[4].upper()
                    info = fields[7]
                    
                    # Combine alternate alleles for variants with the same position
                    if position in variants:
                        variants[position]['alt'].append(alt)
                    else:
                        variants[position] = {'chromosome': chromosome, 'ref': ref, 'alt': [alt], 'info': info}

    # Process the selected variants as per your needs
    # Here, we will print the combined variant information and the IUPAC code if available
    with open(output_file, 'w') as outfile:
        outfile.write("Chromosome\tPosition\tReference\tCombined_Alternate\tCombined_Alleles\tIUPAC\n")

        for position, variant in variants.items():
            chromosome = variant['chromosome']
            ref = variant['ref']
            combined_alt = ','.join(variant['alt'])
            info = variant['info']
            combined_alleles = ref + ',' + combined_alt
            combined_alleles_removedComma = combined_alleles.replace(",", "")
            iupac_code = find_iupac_code(combined_alleles_removedComma)
        
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome, position, ref, combined_alt, combined_alleles, iupac_code))




if __name__ == '__main__':
    args = parse_arguments()
    vcf_file = args.vcf
    output_file = args.out
    parse_vcf(vcf_file, output_file)
