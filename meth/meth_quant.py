import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from natsort import natsorted

data_ih = pd.read_csv("/Users/jon/Documents/data/hapmap/meth/NA12878/IH78cut.1.bedmethyl",
                   names = ['chrom', 'pos', 'code', 'pct'], sep="\t")

data_ont = pd.read_csv("/Users/jon/Documents/data/hapmap/meth/NA12878/ONT78cut.1.bedmethyl",
                   names = ['chrom', 'pos', 'code', 'pct'], sep="\t") 

#Filter rows based on threshold
def meth_filter(data, thresh):
    df = data.drop(data[data.pct < thresh].index)
    return df

#Print h/m/a stats - how many h/m/a;
def hma_stats(data):
    h_subset = data[data['code'].isin(['h'])]
    m_subset = data[data['code'].isin(['m'])]
    a_subset = data[data['code'].isin(['a'])]
    subs_sizes = (h_subset.shape[0], m_subset.shape[0], a_subset.shape[0])
    
    print(f"h : {round(subs_sizes[0] / data.shape[0] * 100, 2)} \
          m : {round(subs_sizes[1] / data.shape[0] * 100, 2)} \
              a : {round(subs_sizes[2] / data.shape[0] * 100, 2)}")
    
hma_stats(meth_filter(data_ih, 20))
hma_stats(meth_filter(data_ont, 20))

#Compare based on position (for how many positions in y, does a call exist in x)
#Loop through data_ih and see if a row with the same chrom and pos exists in data_ont
def call_check(df1, df2, mode=0):
    # Merge on both 'chrom' and 'pos' columns (inner join keeps only matches)
    if mode == 0:
        merged = df1.merge(df2[['chrom', 'pos']], on=['chrom', 'pos'], how='inner')
    #Compare further (how many of these positions that exist, is it the same pct)
    elif mode == 1:
        merged = df1.merge(df2[['chrom', 'pos', 'pct']], on=['chrom', 'pos', 'pct'], how='inner')
    # Drop duplicates in case df2 has multiple matches for the same key
    return merged.drop_duplicates()