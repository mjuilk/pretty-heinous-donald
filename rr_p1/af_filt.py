#!/usr/bin/env python3

#########################
# script for filtering  #
# allele frequency from #
# tab-sep list of vars  #
#########################

import sys
import numpy as np
import pandas as pd
from gnomad_db.database import gnomAD_DB

#establish db connection
loc = "/dss/work/doad5844/gnomad/sqlite3"
db = gnomAD_DB(loc, gnomad_version="v4")

#args for input file and thresh
fin = sys.argv[1]
print(f"Input variants file : {fin}")
thresh = sys.argv[2]

#filtering variants and by AF
with open(fin, 'r') as all_vars:
    all_vars_lns = [ln.split() for ln in all_vars][1:]
    fout = open('af_filt.tsv', 'w')
    fout.write('CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tRSID\tMAF\n')
    print("looping through variants now...")
    for ln in all_vars_lns:
        string = f"{ln[0][3:]}:{ln[1]}:{ln[2]}>{ln[3]}"
        af = db.get_info_from_str(string, "AF")
        joined = '\t'.join(ln)
        if isinstance(af, float):
            if af < float(thresh):
                fout.write(f"{joined}\t{af}\n")
                print(f"##### \t {string} [THRESH MET] : \t {af} \t ##### \n")
            else:
                pass
        else:
            if "," in string:
                print(f"##### \t MULTIALLELIC SITE : \t {string} \t #####")
                sites = string.split(":")[-1].split(">")[-1].split(",")
                for site in sites:
                    string = f"{ln[0][3:]}:{ln[1]}:{ln[2]}>{site}"
                    af = db.get_info_from_str(string, "AF")
                    try:
                        if af < float(thresh):
                            ln[4] = site
                            joined = '\t'.join(ln)
                            print(f"##### \t {string} [THRESH MET] : \t {af} \t ##### \n")
                            fout.write(f"{joined}\t{af}\n")
                        else:
                            print(f"##### \t {string} [THRESH !MET] : \t {af} \t #####")
                    except:
                        print(f"No gnomAD match found! for {string}")
                        fout.write(f"{joined}\n")
    fout.close()
