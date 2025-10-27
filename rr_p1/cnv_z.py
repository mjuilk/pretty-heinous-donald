import numpy as np
import pandas as pd

# Read the TSV file
df = pd.read_csv('/Users/jon/Documents/data/rer/out/c4/test.bed', sep='\t',  dtype={'coverage': float})

# Calculate z-scores for the 'coverage' column
mean = df['coverage'].mean()
stddev = df['coverage'].std()
df['zscore'] = (df['coverage'] - mean) / stddev
thresh = -4

# Filter rows where zscore <= thresh
filtered_df = df[df['zscore'] <= thresh]

# Save the result to a new TSV file
df.to_csv('/Users/jon/Documents/data/rer/out/c4/c4_z.bed', sep='\t', index=False, header=True)
filtered_df.to_csv('/Users/jon/Documents/data/rer/out/c4/c4_z4.bed', sep='\t', index=False, header=True)
