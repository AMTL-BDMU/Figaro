#!/usr/bin/env python3


import csv
import sys
import os

def combine_csv_files(input_files, output_file):
    # Initialize a dictionary to store drug data
    drug_data = {}

    # Loop through each input file
    for input_file in input_files:
        with open(input_file, 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                sample = row['Sample']
                subtype = row['Subtype']
                # Loop through the drug columns in the specified order
                drugs = ['ATV', 'DRV', 'FPV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'ABC', 'AZT', 'D4T', 'DDI', 'FTC', 'LMV', 'TDF', 'DOR', 'EFV', 'ETR', 'NVP', 'RPV', 'BIC', 'CAB', 'DTG', 'EVG', 'RAL']
                for drug in drugs:
                    # If the drug data for this sample is missing, set it to "N/A"
                    drug_data.setdefault((sample, subtype), {})[drug] = row.get(drug, 'N/A')

    # Write the combined data to the output file
    with open(output_file, 'w', newline='') as csv_output:
        fieldnames = ['Sample', 'Subtype'] + drugs
        writer = csv.DictWriter(csv_output, fieldnames=fieldnames)
        writer.writeheader()
        for (sample, subtype), drug_values in drug_data.items():
            row_data = {'Sample': sample, 'Subtype': subtype, **drug_values}
            writer.writerow(row_data)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python combine_csv.py output_file.csv input_file1.csv input_file2.csv ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    if os.path.exists(output_file):
        print("Output file already exists. Please choose a different name.")
        sys.exit(1)

    combine_csv_files(input_files, output_file)
    print(f"Combined data saved to {output_file}")
