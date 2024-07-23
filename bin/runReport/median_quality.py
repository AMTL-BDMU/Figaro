#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from Bio import SeqIO
import gzip

def median_quality_per_position(fastq_file):
    if fastq_file.endswith('.gz'):
        handle = gzip.open(fastq_file, "rt")
    else:
        handle = open(fastq_file, "rt")
    
    records = SeqIO.parse(handle, "fastq")
    quality_scores = []

    for record in records:
        seq_len = len(record)
        if len(quality_scores) < seq_len:
            quality_scores.extend([[] for _ in range(seq_len - len(quality_scores))])
        for i, qual in enumerate(record.letter_annotations["phred_quality"]):
            quality_scores[i].append(qual)
    
    handle.close()
    medians = [np.median(qual) for qual in quality_scores]
    return medians

def write_tsv(medians, output_file):
    with open(output_file, "w") as out:
        out.write("Position\tMedian_Quality\n")
        for position, median_quality in enumerate(medians, start=1):
            out.write(f"{position}\t{median_quality}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate median quality per position in a FASTQ file.")
    parser.add_argument("--inFastq", required=True, help="Input FASTQ file (can be gzipped)")
    parser.add_argument("--outSummary", required=True, help="Output TSV file")
    
    args = parser.parse_args()
    
    medians = median_quality_per_position(args.inFastq)
    write_tsv(medians, args.outSummary)
    
    print(f"Median quality scores per position written to {args.outSummary}")

if __name__ == "__main__":
   main()