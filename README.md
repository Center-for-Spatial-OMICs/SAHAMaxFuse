# SAHA MaxFuse

This repository contains the RNA-protein integration of cell type annotations for the SAHA project.

## PreprocessPRT.ipynb

This script takes the protein Anndata object and splits it into individual organs to prepare for integration. Furthermore, it clusters each individual organ protein Anndata.

## PreprocessRNA.ipynb

This script takes the RNA Anndata object and splits it into individual organs to prepare for integration.

## MaxFuse.py

This script performs the cell matching using the MaxFuse package. Before running the script, there are a few things to edit.

In line 19, change the obs_column to the selected cluster obs name on the protein Anndata object that encompasses the number of cell types from the RNA annotations.
```python
obs_column = 'leiden_3'
```

In line 22-23, change how Anndata objects for protein and RNA are named so that they can be read into the script.
```python
adata_PR = ad.read_h5ad(f'PRT_{tissue}.h5ad')
adata_RNA = ad.read_h5ad(f'RNA_{tissue}.h5ad')
```

In line 33, change what the RNA annotation obs column is named.
```python
labels_rna = adata_RNA.obs['Insitutype_Labelled'].to_numpy()
```

In line 39, change what the csv file is called that contains the analogous protein-RNA name matching. 
```python
correspondence = pd.read_csv('./protein_gene_conversion.csv')
```

Run the script like so and it will generate a matching_{tissue}_{obs_column}.csv file.
```bash
python MaxFuse.py "$tissue"
```

## PostMaxFuse.ipynb

This script takes the matching and assigns cell typing to the protein Anndata object under the obs column name 'celltype_mf'. You can check how the integration performed by examining the cell expression dot plot.

![717e6f09-3524-49f9-8fca-bb98b33a02c2](https://github.com/user-attachments/assets/323cab07-872c-43df-9ddc-1c93c1c2951f)

And the UMAP.

![0f35b18b-3eab-46c0-a6e1-ee426919ffc9](https://github.com/user-attachments/assets/6e7612f7-86f7-4cb7-8cc8-58e1286530ab)

## PostMaxFuseGrannular.ipynb

This script is similar to the one above but instead of integrating broad cell types, we integrated granular cell types plotted the results.

![a95a1660-50e7-4926-ae90-cfdcc26f4524](https://github.com/user-attachments/assets/7b88a37a-beef-457d-ab03-798017745ab9)
![b5da53e9-90f7-41ca-8da0-25311ed7b78e](https://github.com/user-attachments/assets/c3fa6feb-8da4-44f7-b63d-f394927b1952)

## BigUMAP.ipynb

This script combines all the individually matched organs and combine together for analysis. We generated the cell expression heatmap from this script.

![28afb169-b3bd-4828-aea1-4e86c13ab717](https://github.com/user-attachments/assets/00c9f1a5-ec91-4157-8743-7a824b9cfc2f)

## PrettyPictures.ipynb

This script generates the spatial mapping of the cell types. First, we figure out what the fovs so we can select regions to plot and compare.

![55109483-d3a7-4f77-814a-3086b17c6b3f](https://github.com/user-attachments/assets/f234784a-0e24-452f-a95e-d2f1c0dac45b)
![59f592c6-de8a-4058-a3ce-ebe175cd3acd](https://github.com/user-attachments/assets/bf70f95c-b70e-4086-b3b9-fb06aed40762)

From this spatial information, we can examine similar regions in protein and RNA and compare how well the integration performed.

![7a6cb83e-0c62-4c70-a28a-ffcf88ed819b](https://github.com/user-attachments/assets/00fccb2a-35a6-46ec-93a9-0e2db06eff80)
![11c9c925-937b-4f13-bbee-40ae5bd394d6](https://github.com/user-attachments/assets/2e42aec7-e917-4020-8b79-c721dc552ef9)





