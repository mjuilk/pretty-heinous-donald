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
    
    # meth node
    G.add_node(meth_id, type="methylation", beta = row['beta'])
    
    # gene node
    G.add_node(gene, type="gene")
    
    # meth-gene edge
    G.add_edge(meth_id, gene, relationship="methylation_repression", weight=row['anti_corr_score'])

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

def extract_copwin_core(G, verbose=False):
    import random
    import time
    from tqdm import tqdm
    
    H = G.copy()
    removed_nodes = []
    all_nodes = list(H.nodes())
    random.shuffle(all_nodes)
    
    if verbose:
        print(f"Initial graph has {len(H)} nodes")
    
    start_time = time.time()

    pbar = tqdm(total=len(H), desc="Pruning nodes", ncols=80)
    while True:
        progress = False
        for node in all_nodes:
            if node not in H:
                continue

            H_candidate = H.copy()
            H_candidate.remove_node(node)

            is_copwin, _ = is_cop_win_graph(H_candidate)
            if is_copwin:
                removed_nodes.append(node)
                H.remove_node(node)
                progress = True
                if verbose:
                    pbar.update(1)
                break  # restart loop after removing a node

        if not progress:
            break

    pbar.close()
    is_copwin_final, dismantling_order = is_cop_win_graph(H)

    print(f"\nFinal Cop-Win graph has {len(H)} nodes after removing {len(removed_nodes)}")
    print(f"Time taken: {time.time() - start_time:.2f}s")

    # Build annotation dataframe
    node_data = []
    for node in G.nodes:
        in_copwin = node in H
        rank = None
        if in_copwin and dismantling_order:
            rank = dismantling_order[::-1].index(node) + 1  # higher = more resilient

        node_data.append({
            "gene": node,
            "in_copwin": in_copwin,
            "dismantling_rank": rank
        })

    df_annotations = pd.DataFrame(node_data)
    return H, removed_nodes, dismantling_order, df_annotations


G_genes = extract_gene_subgraph(G)

is_copwin, dismantle_order = is_cop_win_graph(G_genes)
if is_copwin:
    print("Cop-Win graph confirmed.")
    print("Dismantling order (most critical genes last):")
    for gene in dismantle_order:
        print(gene)
else:
    print("Graph is not Cop-Win. Consider extracting a Cop-Win subgraph or pruning.")

G_copwin, removed_genes, dismantling_order, annotations = extract_copwin_core(G_genes)

# View top-ranked resilient genes
top_genes = annotations[annotations['in_copwin']].sort_values("dismantling_rank")
print(top_genes.head(10))

is_copwin, dismantle_order = is_cop_win_graph(G_copwin)
if is_copwin:
    print("Cop-Win graph confirmed.")
    print("Dismantling order (most critical genes last):")
    for gene in dismantle_order:
        print(gene)
else:
    print("Graph is not Cop-Win. Consider extracting a Cop-Win subgraph or pruning.")
