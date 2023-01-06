import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#%load_ext rpy2.ipython

## seurat information
sample_obs = pd.read_csv(snakemake.input["cellID_obs"])
umap = pd.read_csv(snakemake.input["cell_embeddings"], header = 0, names = ["Cell ID", "UMAP_1", "UMAP_2"])
celltypes = pd.read_csv(snakemake.input["celltypes"], header = 0, names = ["Cell ID", "celltypes"])
cell_clusters = pd.read_csv(snakemake.input["clusters"], header = 0, names = ["Cell ID", "clusters"])

## loom file
for sample in snakemake.params["samples"]:
    exec(sample + " = anndata.read_loom(snakemake.params['datapath'] + '" + sample + "/velocyto/" + sample + ".loom')")
    exec(sample + ".var_names_make_unique()")
    #exec(sample + ".obs.index = [re.sub(r'^(" + sample + "):(\w*x)$', '\1-\2-1', i) for i in " + sample + ".obs.index]")
    exec(sample + ".obs.index = [re.sub(r':', '-', i) for i in " + sample + ".obs.index]")
    exec(sample + ".obs.index = [re.sub(r'x', '-1', i) for i in " + sample + ".obs.index]")
    exec(sample + " = " + sample + "[np.isin(" + sample + ".obs.index,sample_obs['x'])]")
    #exec("cellID_obs_" + sample + " = sample_obs[cellID_obs_" + sample + "[0].str.contains('" + sample + "-')]")
    #exec(sample + " = " + sample + "[np.isin(" + sample + ".obs.index, cellID_obs_" + sample + ")]")

## merge
exec("myadata = " + snakemake.params['samples'][0] + ".concatenate(" + ','.join(snakemake.params['samples'][1:]) + ")")
myadata.obs.index = [re.sub(r'-1-\d+', '-1', i) for i in myadata.obs.index]


## add umap
myadata_index = pd.DataFrame(myadata.obs.index)
print(myadata_index.head())
myadata_index = myadata_index.rename(columns = {0:'Cell ID'})

#umap = umap_cord.rename(columns = {"Unnamed: 0":"Cell ID"})
umap_ordered = myadata_index.merge(umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
umap_ordered.head()
myadata.obsm['X_umap'] = umap_ordered.values

## add celltypes clusters etc.
celltypes_ordered = myadata_index.merge(celltypes, on = "Cell ID")
celltypes_ordered = celltypes_ordered.iloc[:,1:]
myadata.obs["celltypes"] = celltypes_ordered.values
clusters_ordered = myadata_index.merge(cell_clusters, on = "Cell ID")
clusters_ordered = clusters_ordered.iloc[:, 1:]
myadata.obs["clusters"] = clusters_ordered.values

## run RNA velocity
scv.pp.filter_and_normalize(myadata)
scv.pp.moments(myadata)
scv.tl.velocity(myadata, mode = "stochastic")
scv.tl.velocity_graph(myadata)
myadata.write(snakemake.output[0], compression='gzip')

scv.tl.recover_dynamics(myadata)
scv.tl.velocity(myadata, mode='dynamical')
scv.tl.velocity_graph(myadata)
myadata.write(snakemake.output[1], compression='gzip')

#myadata.filename = snakemake.output[0]
#pdf = PdfPages(snakemake.params["velocity_path"] + "ratio_splicedunspliced.pdf")
#scv.pl.proportions(myadata, show = False)
#plt.savefig(f'{snakemake.params["velocity_path"]}/{ratio_splicedunspliced}.png',dpi=1200)
#pdf.savefig(plt.gcf())
#pdf.close()
#pdf = PdfPages(snakemake.params["velocity_path"] + "embedding.pdf")
#scv.pl.velocity_embedding(myadata, basis = 'umap')
#pdf.savefig()
#pdf.close()
#pdf = PdfPages(snakemake.params["velocity_path"] + "embedding_stream.pdf")
#scv.pl.velocity_embedding_stream(myadata, basis = 'umap')
#pdf.savefig()
#pdf.close()



