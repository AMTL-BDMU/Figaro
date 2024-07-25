#!/usr/bin/env python3

import json
import argparse

def update_sample_result(input_json, output_json, samples, result):
    with open(input_json, 'r') as f:
        data = json.load(f)

    for sample in data['samples']:
        if sample['sample_name'] in samples:
            sample['overallResult'] = result

    with open(output_json, 'w') as f:
        json.dump(data, f, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update sample overall result in JSON file.')
    parser.add_argument('--inJSON', type=str, required=True, help='Input JSON file path')
    parser.add_argument('--outJSON', type=str, required=True, help='Output JSON file path')
    parser.add_argument('--result', type=str, required=True, choices=['Passed', 'Failed'], help='New overall result to set')
    parser.add_argument('--sample', type=str, nargs='+', required=True, help='Sample names to update')

    args = parser.parse_args()

    update_sample_result(args.inJSON, args.outJSON, args.sample, args.result)
