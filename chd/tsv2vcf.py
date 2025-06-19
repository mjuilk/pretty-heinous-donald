#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 12:09:26 2025

code for converting the CHD-CANCER tsv files to VCF format files for Gregor

rsID search is a bit complicated
- for monoallelic variants, it's fine but for multiallelic variants, it encounters issues...
- no current solution for this
- example is rs773340979 where the chr1:g.155900396C>A returns no rsID but chr1:g.155900396C>T does


@author: jon
"""

import os
import myvariant

os.chdir("/Users/jon/Documents/data/chdcan/chd_cancer/analysis")

with open("cancer_germline_LGgenes_list_rare_variants.tsv", "r") as fin, open("cglLG_rare.vcf", "w") as fout:
    mv = myvariant.MyVariantInfo()
    for var in fin.readlines()[1:]:
        spl = var.split()
        geno = spl[0].split(':')
        fetch = mv.getvariant(f"{geno[0]}:g.{geno[1]}{geno[2]}>{geno[3]}", assembly = 'hg38')
        if fetch == None:
            rsid = '.'
        else:
            try:
                rsid = fetch['dbsnp']['rsid']
            except KeyError:
                rsid = '*'
        fout.write(f"{geno[0][3:]}\t{geno[1]}\t{rsid}\t{spl[1]}\t{geno[2]}\t{geno[3]}\n")
        