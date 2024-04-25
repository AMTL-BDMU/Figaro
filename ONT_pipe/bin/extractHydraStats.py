#!/usr/bin/env python3

import argparse
import csv

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Transpose data and add sample name.')
parser.add_argument('--in', dest='input_file', required=True, help='Input text file')
parser.add_argument('--sample', dest='sample_name', required=True, help='Sample name')
parser.add_argument('--out', dest='output_file', required=True, help='Output CSV file')
args = parser.parse_args()

# Read data from the input file
data = {}
with open(args.input_file, 'r') as file:
    for line in file:
        key, value = line.strip().split(': ')
        data[key] = value

# Add the "Sample" header and sample name
data["Sample"] = args.sample_name

# Reorder columns to move "Sample" to the first position
# column_order = ['Sample'] + [col for col in data.columns if col != 'Sample']


data_reordered = {
    'Sample' : data['Sample'],
    'Input Size' : data['Input Size'],
    'Number of reads filtered due to length' : data['Number of reads filtered due to length'],
    'Number of reads filtered due to average quality score' : data['Number of reads filtered due to average quality score'],
    'Number of reads filtered due to presence of Ns' : data['Number of reads filtered due to presence of Ns'],
    'Number of reads filtered due to excess coverage' : data['Number of reads filtered due to excess coverage'],
    'Number of reads filtered due to poor mapping' : data['Number of reads filtered due to poor mapping'],
    'Percentage of reads filtered' : data['Percentage of reads filtered']
}

# Write the transposed data_reordered to a CSV file with a comma delimiter
with open(args.output_file, 'w', newline='') as csv_file:
    csv_writer = csv.DictWriter(csv_file, fieldnames=data_reordered.keys(), delimiter=',')
    csv_writer.writeheader()
    csv_writer.writerow(data_reordered)

print(f"Data has been written to {args.output_file}")

