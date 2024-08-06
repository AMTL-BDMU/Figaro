#!/bin/python3

import json
import argparse

def add_subtype_to_json(in_json, subtype, out_json):
    # Read the input JSON file
    with open(in_json, 'r') as file:
        data = json.load(file)
    
    # Add the subtype to the "subtypeText" field
    if isinstance(data, dict):
        data['subtypeText'] = subtype
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, dict):
                item['subtypeText'] = subtype
    
    # Write the updated data to the output JSON file
    with open(out_json, 'w') as file:
        json.dump(data, file, indent=4)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Add subtype to JSON file.')
    parser.add_argument('--inJSON', type=str, required=True, help='Input JSON file')
    parser.add_argument('--subtype', type=str, required=True, help='Subtype to add')
    parser.add_argument('--outJSON', type=str, required=True, help='Output JSON file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Add the subtype to the JSON file
    add_subtype_to_json(args.inJSON, args.subtype, args.outJSON)

if __name__ == '__main__':
    main()

