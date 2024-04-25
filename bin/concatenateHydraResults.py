#!/usr/bin/env python3

import pandas as pd
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description="Concatenate CSV files after specified markers")

# Add argument for input CSV file after --hydra
parser.add_argument("--hydra", required=True, help="Input CSV file to concatenate after --hydra")

# Add argument for input CSV file after --coverage
parser.add_argument("--coverage", required=True, help="Input CSV file to concatenate after --coverage")

# Add argument for the output concatenated file
parser.add_argument("--out", required=True, help="Output concatenated CSV file")

# Parse the command-line arguments
args = parser.parse_args()

# Read the two input CSV files
file1 = args.hydra
file2 = args.coverage

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Concatenate the dataframes while excluding the second "Sample" column
concatenated_df = pd.concat([df1, df2.iloc[:, 1:]], axis=1)

# Save the concatenated dataframe to the specified output CSV file
concatenated_df.to_csv(args.out, index=False)

