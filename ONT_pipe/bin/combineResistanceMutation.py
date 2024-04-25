#!/usr/bin/env python3

import csv
import sys
import os
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description="Combine multiple CSV files into one.")
parser.add_argument("input_files", nargs="+", help="Input CSV files to combine")
parser.add_argument("--out", required=True, help="Output CSV file")

# Parse the command-line arguments
args = parser.parse_args()

# Extract the input file names from the command-line arguments
input_files = args.input_files

# Use the specified output file
output_file = args.out

# Define the desired column order
desired_columns = [
    "Sample",
    "Gene",
    "Category",
    "Surveillance",
    "Wildtype",
    "Position",
    "Mutation",
    "Mutation Frequency",
    "Coverage"
]
						



# Initialize an empty list to store combined data
combined_data = []

# Initialize an empty list to store headers
headers = []

# Loop through each input file
for file_name in input_files:
    with open(file_name, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        
        # Read the first row (header) of each input file
        file_headers = reader.fieldnames
        
        # Append the headers from this file to the list
        headers.extend(file_headers)

        # Read and append the entries from this file
        for row in reader:
            combined_data.append(row)

# Deduplicate headers (in case there are duplicates)
headers = list(set(headers))


# Ensure the headers are in the desired order
ordered_headers = [col for col in desired_columns if col in headers]

# Write the combined data to the output CSV file with the desired column order
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=ordered_headers)
    
    # Write the header row
    writer.writeheader()
    
    # Write the combined entries
    writer.writerows(combined_data)

print(f'Combined data written to {output_file}')
