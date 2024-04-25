#!/usr/bin/env python3

import json
import argparse

def check_aligned_gene_sequences(input_file, sample_prefix):
    try:
        # Load JSON data from the input file
        with open(input_file, 'r') as file:
            data = json.load(file)

        # Check if the 'alignedGeneSequences' array is empty for the first item in the JSON data
        if len(data) > 0 and not data[0]["alignedGeneSequences"]:
            print("The 'alignedGeneSequences' array is empty.")
            # Create an empty JSON file named 'noAlignment.json'
            with open('noAlignment.json', 'w') as output_file:
                json.dump(data, output_file, indent=4)
        else:
            print("The 'alignedGeneSequences' array is not empty.")
            # Create a JSON file with the specified consensus prefix
            consensus_file = f'consensus_{sample_prefix}.json'
            with open(consensus_file, 'w') as output_file:
                json.dump(data, output_file, indent=4)
            print(f"Consensus JSON file '{consensus_file}' created.")

    except FileNotFoundError:
        print(f"Error: The JSON file '{input_file}' does not exist.")
    except (json.JSONDecodeError, KeyError):
        print(f"Error: Invalid JSON data or missing 'alignedGeneSequences' key.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check if the 'alignedGeneSequences' array is empty in a JSON file.")
    parser.add_argument("--json", required=True, help="Path to the input JSON file")
    parser.add_argument("--sample", required=True, help="Sample prefix for the output file")

    args = parser.parse_args()
    json_file = args.json
    sample_prefix = args.sample

    check_aligned_gene_sequences(json_file, sample_prefix)
