#!/usr/bin/env python

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculate z-scores for coverage data and filter by threshold.")
    parser.add_argument("input_file", help="Path to the input BED file")
    parser.add_argument("threshold", type=float, help="Threshold for filtering z-scores (e.g., -4)")
    args = parser.parse_args()

    # Read the TSV file
    df = pd.read_csv(args.input_file, sep='\t', dtype={'coverage': float})

    # Calculate z-scores for the 'coverage' column
    mean = df['coverage'].mean()
    stddev = df['coverage'].std()
    df['zscore'] = (df['coverage'] - mean) / stddev

    # Filter rows where zscore <= threshold
    filtered_df = df[df['zscore'] <= args.threshold]

    # Save the full DataFrame with z-scores
    df.to_csv(f"{args.input_file.rsplit('.', 1)[0]}_z.bed", sep='\t', index=False, header=True)

    # Save the filtered DataFrame, incorporating the threshold into the filename
    filtered_df.to_csv(f"{args.input_file.rsplit('.', 1)[0]}_z{abs(args.threshold)}.bed", sep='\t', index=False, header=True)

if __name__ == "__main__":
    main()
