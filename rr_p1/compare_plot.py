#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting for REREAD Phase 1 Short-Read Comparison

@author: jon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

rr2_sen = pd.read_csv("/Users/jon/Documents/data/rer/out/rr2/sen/filtered_vars.tsv",
                   sep='\t')

rr2_sarek = pd.read_csv("/Users/jon/Documents/data/rer/out/rr2/sarek//filtvars_merged.tsv",
                   sep='\t')

rr12_sen = pd.read_csv("/Users/jon/Documents/data/rer/out/rr12/sen/filtered_vars.tsv",
                   sep='\t')

rr12_sarek = pd.read_csv("/Users/jon/Documents/data/rer/out/rr12/sarek/filtered_vars.tsv",
                   sep='\t')

#Consequences Pie
def consq_pie(data):
    counts = data['CONSEQUENCE'].value_counts()[:7]
    df = pd.DataFrame({'labels' : list(counts.index), 'values' : list(counts)})
    
    #colors
    colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(6,6))
    plt.pie(df['values'], autopct='%1.1f%%', pctdistance = 1.1, colors = colors)
    #fig.legend(labels=df['labels'])
    #sns.move_legend(fig, "upper left", bbox_to_anchor=(0.8, 1))
    plt.title('Consequence Distribution')
    plt.tight_layout()
    plt.show()
    

#bar plot showing which genes have most variants (use rna rr3 sandbox as template)
def gene_bar(data):
    counts = data['GENE'].value_counts()[:50].sort_index()
    df = pd.DataFrame({'gene' : list(counts.index), 'count' : list(counts)})
    fig, ax = plt.subplots(figsize=(11, 5)) 
    sns.barplot(df, x = 'gene', y = 'count',
                hue = 'count', palette = "flare")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.yticks(np.arange(0,2250, 250))
    #ax.set_title("Variant Distribution Across Genes (50 Highest)")
    ax.set_xlabel("Gene name")
    ax.set_ylabel("Variant count")
    plt.tight_layout()
    plt.show()

#something with visualising rsid distribution?
#maybe for which gene have a bar next to it representing how many NA's occur in 
#the RSID column