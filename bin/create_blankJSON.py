#!/usr/bin/env python3

import json
import argparse
import subprocess
def create_report(title, out_json):
    # Get the current date and time in GMT
    report_date = subprocess.check_output(["date"]).decode().strip()

    # Construct the JSON structure
    report = {
        "report_title": f"Summary Report: {title}",
        "report_date": report_date,
        "samples": []
    }

    # Write to the output JSON file
    with open(out_json, 'w') as json_file:
        json.dump(report, json_file, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a summary report in JSON format.")
    parser.add_argument("--title", required=True, help="Title to append to the report.")
    parser.add_argument("--outJSON", required=True, help="Output JSON file path.")
    
    args = parser.parse_args()
    
    create_report(args.title, args.outJSON)
