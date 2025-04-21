import numpy as np
import pandas as pd
from scipy.io import mmread
import scipy.io

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (6, 4)

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc
import maxfuse as mf

import sys

tissue = sys.argv[1]
print(tissue)
obs_column = 'leiden_1'

# CODE TIME!!
adata_PR = ad.read_h5ad(f'PRT_{tissue}.h5ad')
adata_RNA = ad.read_h5ad(f'RNA_{tissue}.h5ad')

adata_RNA.X = np.asarray(adata_RNA.X.todense())

sc.pp.scale(adata_PR)
sc.pp.scale(adata_RNA)

print(adata_RNA)
print(adata_PR)

labels_rna = adata_RNA.obs['Insitutype_Broad'].to_numpy()
labels_codex = adata_PR.obs[obs_column].to_numpy()

print(labels_rna)
print(labels_codex)

correspondence = pd.read_csv('./protein_gene_conversion.csv')
correspondence.head()

rna_protein_correspondence = []

for i in range(correspondence.shape[0]):
    curr_protein_name, curr_rna_names = correspondence.iloc[i]
    if curr_protein_name not in adata_PR.var_names:
        continue
    if curr_rna_names.find('Ignore') != -1: # some correspondence ignored eg. protein isoform to one gene
        continue
    curr_rna_names = curr_rna_names.split('/') # eg. one protein to multiple genes
    for r in curr_rna_names:
        if r in adata_RNA.var_names:
            rna_protein_correspondence.append([r, curr_protein_name])
            
rna_protein_correspondence = np.array(rna_protein_correspondence)

# Columns rna_shared and protein_shared are matched.
# One may encounter "Variable names are not unique" warning,
# this is fine and is because one RNA may encode multiple proteins and vice versa.
rna_shared = adata_RNA[:, rna_protein_correspondence[:, 0]].copy()
protein_shared = adata_PR[:, rna_protein_correspondence[:, 1]].copy()

rna_shared = rna_shared.X.copy()
protein_shared = protein_shared.X.copy()

# make sure no feature is static
rna_active = adata_RNA.X
protein_active = adata_PR.X
rna_active = rna_active[:, rna_active.std(axis=0) > 1e-5] # these are fine since already using variable features
protein_active = protein_active[:, protein_active.std(axis=0) > 1e-5] # protein are generally variable

# inspect shape of the four matrices
print(rna_active.shape)
print(protein_active.shape)
print(rna_shared.shape)
print(protein_shared.shape)

# call constructor for Fusor object
# which is the main object for running MaxFuse pipeline
fusor = mf.model.Fusor(
    shared_arr1=rna_shared,
    shared_arr2=protein_shared,
    active_arr1=rna_active,
    active_arr2=protein_active,
    labels1=None,
    labels2=None
)

fusor.split_into_batches(
    max_outward_size=8000,
    matching_ratio=4,
    metacell_size=2,
    verbose=True
)

fusor.construct_graphs(
    n_neighbors1=15,
    n_neighbors2=15,
    svd_components1=40,
    svd_components2=15,
    resolution1=2,
    resolution2=2,
    # if two resolutions differ less than resolution_tol
    # then we do not distinguish between then
    resolution_tol=0.1,
    verbose=True
)

fusor.find_initial_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=25, svd_components2=20
)

fusor.refine_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=40, svd_components2=None,
    cca_components=25,
    n_iters=1,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)

fusor.filter_bad_matches(target='pivot', filter_prop=0.5)

pivot_matching = fusor.get_matching(order=(2, 1),target='pivot')

lv1_acc = mf.metrics.get_matching_acc(matching=pivot_matching, 
    labels1=labels_rna, 
    labels2=labels_codex,
    order = (2,1)
)
lv1_acc

# We can inspect the first pivot pair.
[pivot_matching[0][0], pivot_matching[1][0], pivot_matching[2][0]]

cm = confusion_matrix(labels_rna[pivot_matching[0]], labels_codex[pivot_matching[1]])

fusor.propagate(
    svd_components1=40, 
    svd_components2=None, 
    wt1=0.7,
    wt2=0.7,
)

fusor.filter_bad_matches(
    target='propagated',
    filter_prop=0.3
)

full_matching = fusor.get_matching(order=(2, 1), target='full_data')

df = pd.DataFrame(list(zip(full_matching[0], full_matching[1], full_matching[2])), 
            columns = ['mod1_indx', 'mod2_indx', 'score'])
# columns: cell idx in mod1, cell idx in mod2, and matching scores

df.to_csv(f'matching_{tissue}_{obs_column}.csv', index=False)