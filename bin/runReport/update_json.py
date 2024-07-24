#!/usr/bin/env python3

import json
import argparse
import re

# Define valid features and their types
valid_features = {
    'readsNumber_initial': float,
    'readsProportion_passed': float,
    'medianReadLength_initial': float,
    'medianReadLength_final': float,
    'medianDepth': float,
    'genomeCoverage': float,
    'overallResult': str,
    'qc_bp': list,
    'qc_phred': list,
    'alignment_bp': list,
    'alignment_depth': list
}

def parse_value(value, expected_type, feature):
    if expected_type == list:
        # Handle potential range in the list
        if '-' in value:
            # Replace ranges with expanded lists
            range_pattern = re.compile(r'(\d+)-(\d+)')
            matches = range_pattern.findall(value)
            for match in matches:
                start, end = map(int, match)
                value = value.replace(f'{start}-{end}', str(list(range(start, end+1)))[1:-1])
        try:
            # Convert string representation to actual list
            value_list = json.loads(f'[{value}]')
            if not isinstance(value_list, list):
                raise ValueError

            # Ensure the list items are of the correct type
            if feature in ['qc_bp', 'alignment_bp']:
                value_list = [str(item) for item in value_list]
            else:
                value_list = [float(item) if isinstance(item, int) else item for item in value_list]

            return value_list
        except (ValueError, json.JSONDecodeError):
            raise ValueError(f"The value for the list should be a valid list.")
    else:
        parsed_value = expected_type(value)
        if expected_type == float:
            parsed_value = round(parsed_value, 2)
        return parsed_value

def update_json_file(input_json, output_json, sample, feature, value):
    try:
        # Load the JSON data
        with open(input_json, 'r') as file:
            data = json.load(file)
        
        # Check if the 'samples' key is present in the JSON data
        if 'samples' not in data or not isinstance(data['samples'], list):
            raise KeyError("The input JSON does not contain a valid 'samples' list.")
        
        # Find the sample
        sample_found = False
        for sample_entry in data['samples']:
            if sample_entry.get('sample_name') == sample:
                sample_found = True
                # Update or add the feature
                if feature in sample_entry:
                    print(f"Updating {feature} for sample {sample}.")
                else:
                    print(f"Adding {feature} for sample {sample}.")
                sample_entry[feature] = value
                break

        # If the sample is not found, create a new entry
        if not sample_found:
            print(f"Sample {sample} not found. Creating a new entry.")
            new_sample = {'sample_name': sample, feature: value}
            data['samples'].append(new_sample)
        
        # Save the updated JSON data
        with open(output_json, 'w') as file:
            json.dump(data, file, indent=4)
        print("JSON file updated successfully.")
    
    except FileNotFoundError:
        print(f"Error: The file {input_json} was not found.")
    except json.JSONDecodeError:
        print(f"Error: The file {input_json} is not a valid JSON file.")
    except KeyError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Setup argument parser
    parser = argparse.ArgumentParser(description='Update JSON file with sample data.')
    parser.add_argument('--json', required=True, help='Input JSON file path.')
    parser.add_argument('--out', required=True, help='Output JSON file path.')
    parser.add_argument('--sample', required=True, help='Sample name to update or add.')
    parser.add_argument('--feature', required=True, choices=valid_features.keys(), help='Feature to update or add.')
    parser.add_argument('--value', required=True, help='Value to set for the specified feature.')

    # Parse arguments
    args = parser.parse_args()

    # Validate and convert value to appropriate type
    expected_type = valid_features[args.feature]
    try:
        args.value = parse_value(args.value, expected_type, args.feature)
    except ValueError as e:
        parser.error(str(e))

    # Update the JSON file
    update_json_file(args.json, args.out, args.sample, args.feature, args.value)
