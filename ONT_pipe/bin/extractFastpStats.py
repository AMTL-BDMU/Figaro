#!/usr/bin/env python3

import json
import argparse
import pandas as pd

# Create a command-line argument parser
parser = argparse.ArgumentParser(description='Extract information from a JSON file and save it as a CSV.')
parser.add_argument('--json', required=True, help='Path to the JSON file')
parser.add_argument('--output', required=True, help='Path to the CSV output file')
parser.add_argument('--sample', required=True, help='Sample name')

# Parse the command-line arguments
args = parser.parse_args()

# Read the JSON data from the specified file
with open(args.json, 'r') as json_file:
    json_data = json_file.read()

# Parse the JSON data
data = json.loads(json_data)

# Extract information
before_filtering = data["summary"]["before_filtering"]
after_filtering = data["summary"]["after_filtering"]

# Create a DataFrame
df = pd.DataFrame({
    'Sample': [args.sample],
    'BeforeFiltering_TotalReads': [before_filtering['total_reads']],
    'AfterFiltering_TotalReads': [after_filtering['total_reads']],
    'BeforeFiltering_TotalBases': [before_filtering['total_bases']],
    'AfterFiltering_TotalBases': [after_filtering['total_bases']],
    'BeforeFiltering_Q20Bases': [before_filtering['q20_bases']],
    'AfterFiltering_Q20Bases': [after_filtering['q20_bases']],
    'BeforeFiltering_Q30Bases': [before_filtering['q30_bases']],
    'AfterFiltering_Q30Bases': [after_filtering['q30_bases']],
    'BeforeFiltering_Q20Rate': [before_filtering['q20_rate']],
    'AfterFiltering_Q20Rate': [after_filtering['q20_rate']],
    'BeforeFiltering_Q30Rate': [before_filtering['q30_rate']],
    'AfterFiltering_Q30Rate': [after_filtering['q30_rate']],
    'BeforeFiltering_Read1MeanLength': [before_filtering['read1_mean_length']],
    'AfterFiltering_Read1MeanLength': [after_filtering['read1_mean_length']],
    'BeforeFiltering_Read2MeanLength': [before_filtering['read2_mean_length']],
    'AfterFiltering_Read2MeanLength': [after_filtering['read2_mean_length']],
    'BeforeFiltering_GCcontent': [before_filtering['gc_content']],
    'AfterFiltering_GCcontent': [after_filtering['gc_content']]
})

# Save the DataFrame as a CSV file
df.to_csv(args.output, index=False)

print(f"Data saved to {args.output}")

