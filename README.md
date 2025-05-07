# SAHA MaxFuse

This repository contains the RNA-protein integration of cell type annotations for the SAHA project.

## PreprocessPRT.ipynb

This script takes the protein Anndata object and splits it into individual organs to prepare for integration. Furthermore, it clusters each individual organ protein Anndata.

## PreprocessRNA.ipynb

This script takes the RNA Anndata object and splits it into individual organs to prepare for integration.

## MaxFuse.py

This script performs the cell matching using the MaxFuse package. 

```python
obs_column = 'leiden_3'
```
