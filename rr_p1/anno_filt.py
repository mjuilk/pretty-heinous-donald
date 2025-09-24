#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 16:01:57 2025

@author: jon
"""

import numpy as np
import pandas as pd

file_path = "/Users/jon/Documents/data/rer/out/rr87/hap/hq_anno.hg38_multianno.csv"
df = pd.read_csv(file_path)

intervar_filt = df[df['InterVar_automated'] == "Uncertain significance"]
intervar_filt.to_csv('/Users/jon/Documents/data/rer/out/rr87/deep/vus', encoding='utf-8')

df['ClinPred_score']=df.ClinPred_score.replace('.',np.nan).astype(float)
df['ClinPred_score'] = pd.to_numeric(df['ClinPred_score'])
clinpred_filt = df[df['ClinPred_score'] >= 0.5]

df['SIFT_score']=df.SIFT_score.replace('.',np.nan).astype(float)
df['SIFT_score'] = pd.to_numeric(df['SIFT_score'])
sift_filt = df[df['SIFT_score'] <= 0.25]

sift_filt['REVEL']=sift_filt.REVEL.replace('.',np.nan).astype(float)
sift_filt['REVEL'] = pd.to_numeric(sift_filt['REVEL'])
revel_filt = sift_filt[sift_filt['REVEL'] >= 0.25]