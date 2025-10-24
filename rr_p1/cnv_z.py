import numpy as np
import pandas as pd

# Read the TSV file
df = pd.read_csv('/Users/jon/Documents/data/rer/out/c1/chr17_b.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'coverage'])

# Calculate z-scores for the 'coverage' column
mean = df['coverage'].mean()
stddev = df['coverage'].std()
df['zscore'] = (df['coverage'] - mean) / stddev

# Save the result to a new TSV file
df.to_csv('/Users/jon/Documents/data/rer/out/c1/chr17_z.bed', sep='\t', index=False, header=False)
