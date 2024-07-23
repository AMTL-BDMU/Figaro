import json
import argparse
import sys

# Define valid features with their expected types
valid_features = [
    'readsNumber_initial',
    'readsProportion_passed',
    'medianReadLength_initial',
    'medianReadLength_final',
    'medianDepth',
    'genomeCoverage',
    'overallResult',
    'qc_bp',
    'qc_phred',
    'alignment_bp',
    'alignment_depth'
]

def check_samples(input_json, output_file, output_json):
    with open(input_json, 'r') as file:
        data = json.load(file)
    
    samples = data.get('samples', [])
    
    missing_features_report = []
    
    for sample in samples:
        missing_features = []
        for feature in valid_features:
            if feature not in sample:
                missing_features.append(feature)
                sample[feature] = 0
        
        if missing_features:
            missing_features_report.append({
                'sample_name': sample.get('sample_name', 'Unknown'),
                'missing_features': missing_features
            })
    
    with open(output_file, 'w') as file:
        file.write("Sample Name\tMissing Features\n")
        for report in missing_features_report:
            file.write(f"{report['sample_name']}\t{', '.join(report['missing_features'])}\n")
    
    with open(output_json, 'w') as file:
        json.dump(data, file, indent=4)
    
    return bool(missing_features_report)

def main():
    parser = argparse.ArgumentParser(description='Check samples for missing features.')
    parser.add_argument('--inJSON', required=True, help='Input JSON file')
    parser.add_argument('--outTXT', required=True, help='Output text file containing the missing features')
    parser.add_argument('--outJSON', required=True, help='Output JSON file with missing features filled')
    
    args = parser.parse_args()
    
    has_missing_features = check_samples(args.inJSON, args.outTXT, args.outJSON)
    
    if has_missing_features:
        sys.exit(1)

if __name__ == "__main__":
    main()
