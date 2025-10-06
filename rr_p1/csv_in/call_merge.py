#!/usr/bin/env python3

 ################################################
#  Script to merge output from diff callers      #
# 1. define paths for all inputs                 #
# 2. read in data and concat all DFs             #
# 3. remove duplicates to create base df         #
# 4. loop and add caller presence info           #
 ################################################

import os
import pandas as pd

#Define paths
sample = "/user/doad5844/results/reread_p1/rr1"
os.chdir(sample)
deep_path = "rr1_hqde.csv"
hap_path = "rr1_hqha.csv"
sen_path = "rr1_hqse.csv"

#Read CSVs
deep_df = pd.read_csv(deep_path)
print(f"[DEEP]\tNumber of variants : {deep_df.shape[0]}")
hap_df = pd.read_csv(hap_path)
print(f"[HAP]\tNumber of variants : {hap_df.shape[0]}")
sen_df = pd.read_csv(sen_path)
print(f"[SEN]\tNumber of variants : {sen_df.shape[0]}")

#Create dedup base DF
dedup = pd.concat([deep_df, hap_df, sen_df]).drop_duplicates(["Chr", "Start", "Ref", "Alt"])
print(f"[DEDUP]\tNumber of variants : {dedup.shape[0]}")

#Create multi-index for comparison
def make_key(df):
    return df[['Chr', 'Start', 'Ref', 'Alt']].apply(tuple, axis=1)

print("Generating keys")
key_unique = make_key(dedup)
key_dv = make_key(deep_df)
key_hc = make_key(hap_df)
key_st = make_key(sen_df)

#Add indicator columns
print("Adding indicator columns...")
dedup['in_deepvariant'] = key_unique.isin(key_dv).replace({True: 'X', False: ''})
dedup['in_haplotypecaller'] = key_unique.isin(key_hc).replace({True: 'X', False: ''})
dedup['in_sentieon'] = key_unique.isin(key_st).replace({True: 'X', False: ''})

#Write output
print("Writing output with semicolons as delimiter...")
dedup.to_csv("rr1_hqme.csv", sep = ';')
