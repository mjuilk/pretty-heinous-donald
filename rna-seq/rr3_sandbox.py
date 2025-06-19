#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:53:45 2025

RNA sandbox to develop and small useful tools and scripts for working with
RNA data from the REREAD project

data has been produced by the nf-core/rnaseq pipeline -> /star_salmon

@author: jon
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv("/Users/jon/Documents/data/rer/out/rr3/rnaseq/star_salmon/salmon.merged.gene_tpm.tsv",
                   sep='\t')

can_genes = pd.read_csv("/Users/jon/Documents/mats/rsrcs/cancerGeneTru.csv")
can_genes_filt = can_genes[['Hugo Symbol', 'chrom', 'start', 'end',
                            'Oncogene', 'Tumor Supp']]
can_genes_filt = can_genes_filt.dropna()

onc_list = can_genes_filt[can_genes_filt['Oncogene'].isin(['Yes'])]
tsg_list = can_genes_filt[can_genes_filt['Tumor Supp'].isin(['Yes'])]

onc_data = data[data['gene_name'].isin(list(onc_list['Hugo Symbol']))]
tsg_data = data[data['gene_name'].isin(list(tsg_list['Hugo Symbol']))]

fig, ax = plt.subplots(figsize=(11, 5)) 
sns.barplot(onc_data, x = 'gene_name', y = 'REREAD3',
            hue = 'REREAD3', palette = "flare", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_title("Oncogene RNA Sequencing Results")
ax.set_xlabel("Gene name")
ax.set_ylabel("Gene count in TPM")
plt.show()

fig, ax = plt.subplots(figsize=(11, 5)) 
sns.barplot(tsg_data, x = 'gene_name', y = 'REREAD3',
            hue = 'REREAD3', palette = "flare", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_title("Tumour suppressor gene RNA Sequencing Results")
ax.set_xlabel("Gene name")
ax.set_ylabel("Gene count in TPM")
plt.show()
