#!/usr/bin/env python3

 ################################################
#  Script to merge output from diff callers      #
# 1. define paths for all inputs                 #
# 2. read in data and concat all DFs             #
# 3. remove duplicates to create base df         #
# 4. create keys based on pos, ref, alt          #
# 5. add indicator columns for each caller       #
 ################################################

import os
import pandas as pd

def merger_sub(data, sam, mode):
    '''
    [INPUT]: data (self explanatory), sam (sample name) and  mode (either pure or hq) for writing output file name
    to be run inside of the main merger function
    '''

    if mode in ["pu", "hq"]:
        dedup = pd.concat([data['deep'], data['hap'], data['sen']]).drop_duplicates(["Chr", "Start", "Ref", "Alt"])

        print(f"[DEDUP]\tNumber of variants : {dedup.shape[0]}")

        def make_key(df):
            return df[['Chr', 'Start', 'Ref', 'Alt']].apply(tuple, axis=1)

        print("Generating keys")
        key_unique = make_key(dedup)
        key_dv = make_key(data['deep'])
        key_hc = make_key(data['hap'])
        key_st = make_key(data['sen'])

#Add indicator columns
        print("Adding indicator columns...")
        dedup['in_deepvariant'] = key_unique.isin(key_dv).replace({True: 'X', False: ''})
        dedup['in_haplotypecaller'] = key_unique.isin(key_hc).replace({True: 'X', False: ''})
        dedup['in_sentieon'] = key_unique.isin(key_st).replace({True: 'X', False: ''})

#Filtering benign and likely benign before writing to csv
        intervar_filt = dedup[~((dedup['InterVar_automated'] == ".") | (dedup['InterVar_automated'] == "Benign") | (dedup['InterVar_automated'] == "Likely benign"))]
        clnsg_filt = multifilt = dedup[((dedup['CLNSIG'] == "Uncertain_significance") | (dedup['CLNSIG'] == "Pathogenic") | (dedup['CLNSIG'] == "other"))]
        cat_filt = pd.concat([intervar_filt, clnsg_filt]).drop_duplicates(["Chr", "Start", "Ref", "Alt"])

#Write output
        print("Writing output with semicolons as delimiter...")
        if mode == "pu":
            out = sam + "_pume.csv"
        elif mode == "hq":
            out = sam + "_hqme.csv"
        cat_filt.to_csv(out, sep = ';')
    else:
        print("NOT A VALID MODE")
        pass

def merger(sam):
    '''
    [INPUT]: 'sam' (str) relative directory and also acts as the sample name for filename insertion
    change directory to the sample directory
    create two dictionaries with the csvs from annovar, one for pure and one for hq
    create the base df (concatenated and deduplicated)
    create keys based on Chr, Start, Ref, and Alt and then add indicator columns
    write output csv
    '''

    os.chdir(sam)
    pu_dict = {
            'deep' : pd.read_csv(sam + "_pude.csv"),
            'hap' : pd.read_csv(sam + "_puha.csv"),
            'sen' : pd.read_csv(sam + "_puse.csv")
            }
    print("=====PURE/NO QC=====\n")
    print(f"[DEEP]\tNumber of variants : {pu_dict['deep'].shape[0]}")
    print(f"[HAP]\tNumber of variants : {pu_dict['hap'].shape[0]}")
    print(f"[SEN]\tNumber of variants : {pu_dict['sen'].shape[0]}")

    merger_sub(pu_dict, sam, "pu")

    hq_dict = {
            'deep' : pd.read_csv(sam + "_hqde.csv"),
            'hap' : pd.read_csv(sam + "_hqha.csv"),
            'sen' : pd.read_csv(sam + "_hqse.csv")
            }

    print("=====HIGH QC=====\n")
    print(f"[DEEP]\tNumber of variants : {hq_dict['deep'].shape[0]}")
    print(f"[HAP]\tNumber of variants : {hq_dict['hap'].shape[0]}")
    print(f"[SEN]\tNumber of variants : {hq_dict['sen'].shape[0]}")

    merger_sub(hq_dict, sam, "hq")

samples = [f"rr{i}" for i in range(1,22)]

for sam in samples:
    os.chdir("/user/doad5844/results/reread_p1")
    print(sam)
    merger(sam)
