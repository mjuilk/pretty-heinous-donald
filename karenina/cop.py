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
import numpy as np

# Initialise graph

G = nx.Graph()

# Load all input data

## Variants

variants_df = pd.read_csv("~/Documents/data/karenina/afm_cancer.tsv", sep="\t",
                          names = ['chr', 'pos', 'ref', 'alt', 'gene', 'csq', 'rsid', 'maf'])
variants_df = variants_df.iloc[1:]
variants_df['maf'] = pd.to_numeric(variants_df['maf'], errors='coerce')
epsilon = 1e-6
variants_df['inv_maf'] = 1 / (variants_df['maf'] + epsilon)

# Min-max normalization to [0, 1] range
min_val = variants_df['inv_maf'].min()
max_val = variants_df['inv_maf'].max()
variants_df['weight'] = (variants_df['inv_maf'] - min_val) / (max_val - min_val)

# Clean up
variants_df.drop(columns=['inv_maf'], inplace=True)
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


# Calculate anti-correlation and repression

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
merged['beta'] = pd.to_numeric(merged['beta']) / 100
merged['norm_expr'] = pd.to_numeric(merged['norm_expr']) * 10

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

for _, row in variants_df.iterrows():
    variant_id = f"{row['rsid']}_{row['chr']}:{row['pos']}"
    gene = row['gene']
    
    # variant node
    G.add_node(variant_id, type="variant", consequence = row['csq'])
    
    # gene node
    G.add_node(gene, type="gene")
    
    # Add variant-gene edge
    G.add_edge(variant_id, gene, relationship="variant_in_gene", weight = row['weight'])
    
for _, row in merged.iterrows():
    meth_id = f"cg_{row['Chromosome']}:{row['Start']}"
    gene = row['gene']
    
    if row['beta'] > 0.6 and row['norm_expr'] < 0.4:
        row['repr_score'] = (row['beta'] - 0.5) * (0.5 - row['norm_expr']) * 4  # emphasis on repression
    else:
        row['repr_score'] = 0
    
    # meth node
    G.add_node(meth_id, type="methylation", beta = row['beta'])
    
    # gene node
    G.add_node(gene, type="gene")
    
    # meth-gene edge
    G.add_edge(meth_id, gene, relationship="methylation_repression", weight=row['repr_score'])

genes_in_graph = [n for n in G.nodes if G.nodes[n].get('type') == 'gene']
ppi_dict = get_string_interactions(genes_in_graph)

for gene, partners in ppi_dict.items():
    for partner, score in partners:
        # Ensure both gene nodes are added with correct type
        for g in [gene, partner]:
            if g not in G:
                G.add_node(g, type="gene")
            elif G.nodes[g].get("type") is None:
                G.nodes[g]["type"] = "gene"

        if G.has_edge(gene, partner):
            G[gene][partner]['score'] = max(G[gene][partner].get('score', 0), score)
        else:
            G.add_edge(gene, partner, score=score)

def get_graph_methyl_score(G, gene):
    score = 0
    for neighbor in G.neighbors(gene):
        if G.nodes[neighbor].get("type") == "methylation":
            edge_data = G[gene][neighbor]
            score += edge_data.get("weight", 0)
    return score

def get_graph_variant_score(G, gene):
    score = 0
    for neighbor in G.neighbors(gene):
        if G.nodes[neighbor].get("type") == "variant":
            edge_data = G[gene][neighbor]
            score += edge_data.get("weight", 0)
    return score


def dismantle_full_graph(G):
    """
    Greedy dismantling adapted for heterogeneous graph.
    Outputs: dismantling order and gene-specific removal ranks.
    """
    dismantling_order = []
    H = G.copy()
    removal_rank = {}
    removal_step = 0

    while H.nodes:
        found = False
        for node in list(H.nodes):
            neighbors = set(H.neighbors(node))
            for other in H.nodes:
                if other == node:
                    continue
                if neighbors.issubset(set(H.neighbors(other)).union({other})):
                    dismantling_order.append(node)
                    if G.nodes[node].get("type") == "gene":
                        removal_rank[node] = removal_step
                        removal_step += 1
                    H.remove_node(node)
                    found = True
                    break
            if found:
                break
        if not found:
            break  # Can't continue dismantling

    return dismantling_order, removal_rank

def compute_gene_scores_graph(G, removal_rank):
    scores = {}
    max_rank = max(removal_rank.values()) if removal_rank else 1

    for node in G.nodes:
        if G.nodes[node].get("type") != "gene":
            continue

        rank_score = 1 - (removal_rank.get(node, max_rank) / max_rank)

        var_score = get_graph_variant_score(G, node)
        meth_score = get_graph_methyl_score(G, node)

        # Normalize
        var_norm = np.log1p(var_score) / 10
        meth_norm = meth_score / 10

        total = rank_score + var_norm + meth_norm

        scores[node] = {
            "dismantling_rank": rank_score,
            "variant_score": var_norm,
            "methylation_score": meth_norm,
            "total_score": total
        }

    return scores

dismantling_order, removal_rank = dismantle_full_graph(G)
gene_scores = compute_gene_scores_graph(G, removal_rank)

# Top genes by total score
top_genes = sorted(gene_scores.items(), key=lambda x: -x[1]['total_score'])
for gene, scores in top_genes[:20]:
    print(f"{gene}: {scores} \n")

