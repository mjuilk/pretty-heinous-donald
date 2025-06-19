import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from natsort import natsorted

data_ih = pd.read_csv("/Users/jon/Documents/data/hapmap/meth/NA12878/IH78cut.1.bedmethyl",
                   names = ['chrom', 'pos', 'code', 'pct'], sep="\t")

data_ont = pd.read_csv("/Users/jon/Documents/data/hapmap/meth/NA12878/ONT78cut.1.bedmethyl",
                   names = ['chrom', 'pos', 'code', 'pct'], sep="\t") 

data = data_ih
h_subset = data[data['code'].isin(['h'])]
m_subset = data[data['code'].isin(['m'])]

def chr_cat(data):
    chr_ls = natsorted(set(data['chrom']))
    x = data['chrom']
    x = np.array(pd.Categorical(x, categories=chr_ls, ordered=True).codes) 
    x += 1
    x_labels = [str(i) for i in list(set(x))]
    x_labels[-2] = 'X'
    x_labels[-1] = 'Y'
    return (x, x_labels)

#color maps
##sequential : 'plasma_r, 'hot_r'
##diverging : 'coolwarm', 'bwr'

fig, ax = plt.subplots(1,1, figsize=(10, 6))
hb = ax.hexbin(chr_cat(data)[0], data['pos'], C=data['pct'],
                  gridsize=50, reduce_C_function=np.mean, cmap='coolwarm')

# Add colorbar
cb = plt.colorbar(hb, ax = ax)
cb.set_label("Mean mod %")

# Set labels and title
ax.set_xlabel("Chromosome")
ax.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
ax.set_xticklabels(chr_cat(h_subset)[1])
ax.set_yticks(np.linspace(0, max(data['pos']), 10))
ax.set_ylabel("Base Position")
ax.set_title("Methylation Map (NA12878_ONT) - 50 bins")

#Show base plot
plt.show()

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(13, 10))
hb_h = ax1.hexbin(chr_cat(h_subset)[0], h_subset['pos'], C=h_subset['pct'],
                  gridsize=50, reduce_C_function=np.mean, cmap="coolwarm")

# Add colorbar
cb = plt.colorbar(hb_h, ax=ax1)
cb.set_label("Mean mod %")

# Set labels and title
ax1.set_xlabel("Chromosome")
ax1.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
ax1.set_xticklabels(chr_cat(h_subset)[1])
ax1.set_yticks(np.linspace(0, max(h_subset['pos']), 10))
ax1.set_ylabel("Base Position")
ax1.set_title("Methylation h_Map (NA12878_ONT) - 50 bins")

hb_m = ax2.hexbin(chr_cat(m_subset)[0], m_subset['pos'], C=m_subset['pct'],
                  gridsize=50, reduce_C_function=np.mean, cmap="coolwarm")

# Add colorbar
cb = plt.colorbar(hb_m, ax=ax2)
cb.set_label("Mean mod %")

# Set labels and title
ax2.set_xlabel("Chromosome")
ax2.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
ax2.set_xticklabels(chr_cat(m_subset)[1])
ax2.set_yticks(np.linspace(0, max(m_subset['pos']), 10))
ax2.set_ylabel("Base Position")
ax2.set_title("Methylation m_Map (NA12878_ONT) - 50 bins")

# Show h and m subplots
plt.show()
