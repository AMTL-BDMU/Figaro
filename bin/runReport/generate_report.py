#!/usr/bin/env python3

import argparse
import json
from jinja2 import Environment, FileSystemLoader


def generate_report(input_json, template_html, report_html):

    # Load the sample data from JSON file
    with open(input_json, 'r') as f: ###################### data json
        data = json.load(f)

    # Extract data for the bar plot and line plot
    sample_names = [sample['sample_name'] for sample in data['samples']]
    sample_genomeCoverages = [sample['genomeCoverage'] for sample in data['samples']]

    sample_readsNumber_initial = [sample['readsNumber_initial'] for sample in data['samples']]
    sample_readsProportion_passed = [sample['readsProportion_passed'] for sample in data['samples']]
    sample_meanReadLength_initial = [sample['meanReadLength_initial'] for sample in data['samples']]
    sample_meanReadLength_final = [sample['meanReadLength_final'] for sample in data['samples']]
    sample_meanDepth = [sample['meanDepth'] for sample in data['samples']]

    readQC_data = [{
        'name': sample['sample_name'],
        'qc_bp': sample['qc_bp'],
        'qc_phred': sample['qc_phred']
    } for sample in data['samples']]

    alignment_data = [{
        'name': sample['sample_name'],
        'alignment_bp': sample['alignment_bp'],
        'alignment_depth': sample['alignment_depth']
    } for sample in data['samples']]

    # Prepare the data for rendering in the template
    template_data = {
        'report_title': data['report_title'],
        'report_date': data['report_date'],
        'samples': data['samples'],
        'sample_names': json.dumps(sample_names),
        'sample_genomeCoverages': json.dumps(sample_genomeCoverages),
        'readQC_data': json.dumps(readQC_data),
        'alignment_data': json.dumps(alignment_data),
        'sample_readsNumber_initial': json.dumps(sample_readsNumber_initial),
        'sample_readsProportion_passed': json.dumps(sample_readsProportion_passed),
        'sample_meanReadLength_initial': json.dumps(sample_meanReadLength_initial),
        'sample_meanReadLength_final': json.dumps(sample_meanReadLength_final),
        'sample_meanDepth': json.dumps(sample_meanDepth)
    }

    # Render the HTML template with the data
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template(template_html) ############## html template
    output_html = template.render(template_data)

    # Save the rendered HTML to a file
    with open(report_html, 'w') as f: ############## output report html
        f.write(output_html)

    print("HTML report generated successfully.")


if __name__ == "__main__":
    # Setup argument parser
    parser = argparse.ArgumentParser(description='Generate report from the json file.')
    parser.add_argument('--json', required=True, help='Input JSON file path.')
    parser.add_argument('--template', required=True, help='Input HTML Template file path.')
    parser.add_argument('--report', required=True, help='Output HTML Report file path.')

    # Parse arguments
    args = parser.parse_args()

    # Update the JSON file
    generate_report(args.json, args.template, args.report)