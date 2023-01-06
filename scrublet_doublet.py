#!/data/wangj/software/Python-3.6.8-install/bin/python3.6

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

if len(sys.argv) > 3:
    print("only first two paramers are used")

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

input_dir = sys.argv[1]
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))

#print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
#print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
scrub.plot_histogram()
plt.savefig("results/doublet/" + sys.argv[2] + "_doublet_histogram.png")

np.savetxt("results/doublet/" + sys.argv[2] + "_doublet_scores.txt",doublet_scores)
np.savetxt("results/doublet/" + sys.argv[2] + "_predicted_doublets.txt",predicted_doublets)

#print('Running UMAP...')
#scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
#print('Done.')
#scrub.plot_embedding('UMAP', order_points=True)
#plt.savefig("results/doublet/" + sys.argv[2] + "_plot_doublet_clustering.png")



