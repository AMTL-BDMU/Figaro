#!/usr/bin/env python3

import json
import argparse
import os

def combine_json_files(input_files, output_file):
    combined_data = {
        "report_title": "Summary Report",
        "report_date": None,  # This could be set to the current date/time if needed
        "samples": []
    }

    for file in input_files:
        with open(file, 'r') as f:
            data = json.load(f)
            combined_data["samples"].extend(data["samples"])
            if combined_data["report_date"] is None:
                combined_data["report_date"] = data["report_date"]

    # Sort the samples based on the 'sample_name' key
    combined_data["samples"] = sorted(combined_data["samples"], key=lambda x: x["sample_name"])

    with open(output_file, 'w') as f:
        json.dump(combined_data, f, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Combine multiple JSON files into one JSON file.")
    parser.add_argument('-i', '--input', nargs='+', required=True, help="Input JSON files")
    parser.add_argument('-o', '--output', required=True, help="Output JSON file")

    args = parser.parse_args()

    if not args.input:
        print("No input files provided.")
        return

    for file in args.input:
        if not os.path.isfile(file):
            print(f"Input file {file} does not exist.")
            return

    combine_json_files(args.input, args.output)

if __name__ == "__main__":
    main()
