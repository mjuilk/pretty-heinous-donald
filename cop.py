#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 11:23:57 2025

Cop-Win graph and dismantling code for KARENINA

@author: jon
"""

# Packages

import pyranges as pr
import pandas as pd
import requests
import time
import networkx as nx

# Initialise graph

G = nx.Graph()

# Load all input data

## Variants

variants_df = pd.read_csv("~/Documents/data/karenina/afm_cancer.tsv", sep="\t",
                          names = ['chr', 'pos', 'ref', 'alt', 'gene', 'csq', 'rsid', 'maf'])
variants_df = variants_df.iloc[1:]
variants_df['weight'] = pd.to_numeric(variants_df['maf']) * 1000  # Scale for visibility
##########CHANGE THIS SO THAT LOWER MAF MEANS HIGHER WEIGHT!!!!!#############
variants_df = variants_df[~variants_df.alt.str.contains(",")] # Filter multiallelic rows


## TSS

pro_df = pd.read_csv("~/Documents/src/karenina/promote", sep="\t",
                       names = ['chr', 'start', 'end', 'gene'])
pro_df = pro_df[pro_df['gene'].isin(list(variants_df['gene']))]

## Mods

mod_df = pd.read_csv("~/Documents/data/karenina/rr3_chr13_17.bedmethyl", sep="\t",
                     names = ['chr', 'pos', 'code', 'beta'])

## Expression

expression_df = pd.read_csv("~/Documents/data/karenina/expression.tsv", sep="\t",
                            names=["gene", "expr"])

expression_df["norm_expr"] = (expression_df["expr"] - expression_df["expr"].min()) / \
                                  (expression_df["expr"].max() - expression_df["expr"].min())


# Calculate anti-correlation

mod_df['Start'] = mod_df['pos'].astype(int)
mod_df['End'] = mod_df['Start'] + 1  # 1bp interval
mod_df = mod_df.rename(columns={'chr': 'Chromosome'})
mod_df = mod_df[['Chromosome', 'Start', 'End', 'beta']]
mod_pr = pr.PyRanges(mod_df)

pro_df = pro_df.rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'})
pro_df = pro_df[['Chromosome', 'Start', 'End', 'gene']]  # Keep only necessary columns
pro_pr = pr.PyRanges(pro_df)

overlap = mod_pr.join(pro_pr).df
overlap['gene'] = overlap['gene'].astype(str)
expression_df.index = expression_df.index.astype(str)
expression_df = expression_df.set_index('gene')
merged = overlap.merge(expression_df[['norm_expr']], left_on='gene', right_index=True, how='left')
merged['anti_corr_score'] = -4 * (merged['norm_expr'] - 0.5) * (merged['beta'] - 0.5)

# Pull PPI

def get_string_interactions(gene_list, species=9606, required_score=400):
    base_url = "https://string-db.org/api/json/network"
    interactions_dict = {}

    for gene in gene_list:
        params = {
            "identifiers": gene,
            "species": species,
            "required_score": required_score
        }
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            print(f"Failed to fetch data for {gene}")
            interactions_dict[gene] = []
            continue

        interactions = response.json()
        partners = []

        for interaction in interactions:
            # STRING uses protein identifiers, which may differ; filter both ways
            partner = interaction['preferredName_B'] if interaction['preferredName_A'] == gene else interaction['preferredName_A']
            score = interaction['score']
            partners.append((partner, score))

        interactions_dict[gene] = partners

        time.sleep(1)

    return interactions_dict

# Network construction

def build_multiomic_graph(variants_df, meth_expr_df, ppi_dict):
    G = nx.Graph()

    # === Add Gene Nodes ===
    genes = set(variants_df['gene']).union(ppi_dict.keys())
    for gene in genes:
        G.add_node(gene, type='gene')

    # === Add Variant Nodes & Connect to Genes ===
    for _, row in variants_df.iterrows():
        variant_node = f"var:{row['chr']}:{row['pos']}:{row['ref']}>{row['alt']}"
        G.add_node(variant_node, type='variant', consequence=row['consequence'], maf=row['maf'])
        G.add_edge(variant_node, row['gene'], relation='variant_to_gene')

    # === Add Methylation Nodes & Connect to Genes ===
    for idx, row in meth_expr_df.iterrows():
        meth_node = f"meth:{row['chr']}:{row['Start']}-{row['End']}"
        G.add_node(meth_node, type='methylation', anticorrelation=row['anti_corr_score'])
        closest_gene = find_closest_gene(row, variants_df)  # function below
        if closest_gene:
            G.add_edge(meth_node, closest_gene, relation='meth_to_gene', weight=abs(row['anticorrelation']))

    # === Add PPI Edges Between Genes ===
    for gene, partners in ppi_dict.items():
        for partner in partners:
            if gene in G and partner in G:
                G.add_edge(gene, partner, relation='ppi', weight=1)

    return G

def find_closest_gene(meth_row, variants_df):
    # Optional: improve with gene TSS positions
    chrom_genes = variants_df[variants_df['chr'] == meth_row['chr']]
    if chrom_genes.empty:
        return None
    gene_dists = chrom_genes.copy()
    gene_dists['dist'] = abs(gene_dists['pos'] - meth_row['start'])
    closest = gene_dists.sort_values('dist').iloc[0]
    return closest['gene']

def is_cop_win_graph(G):
    # A graph is Cop-Win iff it is dismantlable
    # We'll use a greedy dismantling check
    dismantling_order = []
    H = G.copy()

    while H.nodes:
        found = False
        for node in list(H.nodes):
            neighbors = set(H.neighbors(node))
            for other in H.nodes:
                if other == node:
                    continue
                if neighbors.issubset(set(H.neighbors(other)).union({other})):
                    dismantling_order.append(node)
                    H.remove_node(node)
                    found = True
                    break
            if found:
                break
        if not found:
            return False, []
    return True, dismantling_order[::-1]  # reverse gives from last removed to first

def extract_gene_subgraph(G):
    # Consider only gene nodes and PPI edges
    gene_nodes = [n for n, d in G.nodes(data=True) if d['type'] == 'gene']
    return G.subgraph(gene_nodes).copy()

