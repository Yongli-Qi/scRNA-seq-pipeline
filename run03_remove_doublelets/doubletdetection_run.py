##############################################
#   annotate cells without doublets
##############################################
import sys
sys.path.append('/home/yq238/.local/lib/python3.11/site-packages')
sys.path.append('/vast/palmer/apps/services/ood/mccleary/var_www_ood_apps/conda/ycrc_default/lib/python3.11/site-packages')
import numpy as np
import doubletdetection
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from sys import argv
#sc.settings.n_jobs=16

_, file_in, file_in2, file_out, file_out2 = argv

adata = sc.read_10x_mtx(
    file_in,
    var_names='gene_symbols', 
    cache=False
)

file_path = file_in2
with open(file_path, 'r') as f:
    real_cells = f.read().splitlines()


adata.obs.index = adata.obs.index.astype(str)
adata = adata[adata.obs.index.isin(real_cells)].copy()
sc.pp.filter_genes(adata, min_cells=1)

clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)
doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
doublet_score = clf.doublet_score()

adata.obs["doublet"] = doublets
adata.obs["doublet_score"] = doublet_score


mask = adata.obs['doublet'] == 0.0
filtered_adata = adata[mask, :]
cell_names = filtered_adata.obs.index
with open(file_out, 'w') as file:
    for name in cell_names:
        file.write(name + '\n')

row_names = adata.obs.index
df_output = pd.DataFrame({
'RowName': row_names,
'DoubletScore': doublet_score
})
df_output.to_csv(file_out2, sep='\t', index=False)

with open('log_file.txt', 'a') as log_file:
    log_file.write('Finished: '+file_in+'\n')



