"""
scGeneset/core/loader.py
"""
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from typing import Optional


def load_data(
   annotation_file: Optional[str] = None,
   group_file: Optional[str] = None,
   raw_file: Optional[str] = None,
   adata_file: Optional[str] = None,
   annotation_col: str = 'final_annotation',
   group_col: str = 'group',
   min_genes: Optional[int] = None,
   min_cells: Optional[int] = None,
   target_sum: Optional[float] = None,
   n_top_genes: Optional[int] = None,
   log_normalize: bool = False,
   verbose: bool = True
) -> anndata.AnnData:
   """
   Load and organize data, perform basic preprocessing, and correctly set .raw attribute.

   Parameters:
       annotation_file: Path to cell annotation file (optional)
       group_file: Path to group file (optional)
       raw_file: Path to raw expression matrix file (optional)
       adata_file: Path to an existing h5ad file (preferred, if provided, others are ignored)
       annotation_col: Column name in adata.obs for cell type annotation, default 'final_annotation'
       group_col: Column name in adata.obs for group information, default 'group'
       min_genes: Minimum number of genes expressed for a cell to be kept (for initial filtering)
       min_cells: Minimum number of cells a gene must be expressed in to be kept (for initial filtering)
       target_sum: Target sum for total counts per cell (for normalization, e.g., 1e4)
       n_top_genes: Number of highly variable genes to select
       log_normalize: Whether to apply log1p transformation after normalization
       verbose: Whether to print detailed messages

   Returns:
       adata: AnnData object containing annotation, group information, and preprocessed data.
   """
   if verbose:
       print("============================================")
       print("Step 1: Loading data")
       print("============================================")

   if adata_file is not None:
       adata = sc.read_h5ad(adata_file)
       if verbose:
           print(f"Data loaded from {adata_file}")
           print(f"  Cells: {adata.n_obs}")
           print(f"  Genes: {adata.n_vars}")
           print(f"  obs columns: {adata.obs.shape[1]}")

       if annotation_col not in adata.obs.columns:
           raise ValueError(f"Annotation column '{annotation_col}' not found in adata.obs.")
       if group_col not in adata.obs.columns:
           raise ValueError(f"Group column '{group_col}' not found in adata.obs.")

   else:
       if annotation_file is None or group_file is None or raw_file is None:
           raise ValueError("Must provide either adata_file or (annotation_file, group_file, raw_file)")

       if annotation_file.endswith(('.csv', '.txt', '.tsv')):
           annot = pd.read_csv(annotation_file, index_col=0, sep=',' if annotation_file.endswith('.csv') else '\t')
       else:
           raise ValueError("Annotation file format not supported (must be .csv, .txt, or .tsv)")

       if group_file.endswith(('.csv', '.txt', '.tsv')):
           group = pd.read_csv(group_file, index_col=0, sep=',' if group_file.endswith('.csv') else '\t')
       else:
           raise ValueError("Group file format not supported (must be .csv, .txt, or .tsv)")

       raw_df = None
       if raw_file.endswith(('.csv', '.txt', '.tsv')):
           raw_df = pd.read_csv(raw_file, index_col=0, sep=',' if raw_file.endswith('.csv') else '\t')
       elif raw_file.endswith('.h5ad'):
           raw_adata_temp = sc.read_h5ad(raw_file)
           raw_df = pd.DataFrame(
               raw_adata_temp.X.toarray() if hasattr(raw_adata_temp.X, 'toarray') else raw_adata_temp.X,
               index=raw_adata_temp.obs_names,
               columns=raw_adata_temp.var_names
           )
       else:
           raise ValueError("Expression matrix file format not supported (must be .csv, .txt, .tsv, or .h5ad)")

       raw_df.index = raw_df.index.astype(str)
       raw_df.columns = raw_df.columns.astype(str)
       annot.index = annot.index.astype(str)
       group.index = group.index.astype(str)

       common_cells = annot.index.intersection(group.index).intersection(raw_df.index)
       if len(common_cells) == 0:
           raise ValueError("No common cell IDs found across annotation, group, and expression matrix.")
       if verbose: print(f"Found {len(common_cells)} common cells.")

       adata = anndata.AnnData(X=raw_df.loc[common_cells].values)
       adata.obs_names = common_cells
       adata.var_names = raw_df.columns

       adata.obs[annotation_col] = annot.loc[common_cells].iloc[:, 0].values
       adata.obs[group_col] = group.loc[common_cells].iloc[:, 0].values

   if adata.raw is None:
       adata.raw = adata.copy()
   elif not isinstance(adata.raw, anndata.AnnData):
       if hasattr(adata.raw, 'to_adata'):
           adata.raw = adata.raw.to_adata()
       else:
           print("Warning: adata.raw exists but is not an AnnData object and lacks .to_adata(). Attempting copy of current adata.")
           adata.raw = adata.copy()

   adata.obs['cell_type'] = adata.obs[annotation_col].astype('category')
   adata.obs['group'] = adata.obs[group_col].astype('category')

   if verbose:
       print(f"Data organization complete.")
       print(f"  - Cell types: {adata.obs['cell_type'].nunique()} ({', '.join(adata.obs['cell_type'].cat.categories.tolist())})")
       print(f"  - Groups: {adata.obs['group'].nunique()} ({', '.join(adata.obs['group'].cat.categories.tolist())})")

   if verbose: print("============================================")
   if verbose: print("Step 2: Basic preprocessing")
   if verbose: print("============================================")

   initial_cells = adata.n_obs
   initial_genes = adata.n_vars

   if min_genes is not None:
       sc.pp.filter_cells(adata, min_genes=min_genes)
       if verbose: print(f"  - Filtered cells with < {min_genes} genes. Remaining cells: {adata.n_obs}")

   if min_cells is not None:
       sc.pp.filter_genes(adata, min_cells=min_cells)
       if verbose: print(f"  - Filtered genes with < {min_cells} cells. Remaining genes: {adata.n_vars}")

   if target_sum is not None:
       sc.pp.normalize_total(adata, target_sum=target_sum)
       if verbose: print(f"  - Normalized total counts to {target_sum} per cell.")

   if log_normalize:
       sc.pp.log1p(adata)
       if verbose: print(f"  - Applied log1p transformation.")

   if n_top_genes is not None:
       sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
       adata = adata[:, adata.var['highly_variable']].copy()
       if verbose: print(f"  - Selected {adata.n_vars} highly variable genes.")

   if verbose:
       print(f"Preprocessing complete.")
       print(f"  Final cells: {adata.n_obs} (reduced from {initial_cells})")
       print(f"  Final genes: {adata.n_vars} (reduced from {initial_genes})")
       print(f"  obs columns: {adata.obs.shape[1]}")

   return adata


def load_pathway_genes(pathway_file):
   """
   Load pathway gene sets from a CSV file or predefined pathway name

   Parameters:
       pathway_file : str (REQUIRED)
           Either:
           - A built-in pathway name: 'inflammation', 'PCD', or 'PTM'
           - A file path to custom CSV file: 'path/to/custom_pathway.csv'

   Returns:
       pathway_dict : dict
           Dictionary of {pathway_name: gene_set}

   Examples:
       pathway_dict = scg.load_pathway_genes('inflammation')
       pathway_dict = scg.load_pathway_genes('/path/to/my_pathways.csv')

   Raises:
       ValueError: If pathway_file is not provided
       FileNotFoundError: If specified file doesn't exist
   """

   if pathway_file is None:
       raise ValueError(
           "Error: 'pathway_file' is required!\n"
           "Please specify either:\n"
           "  - Built-in pathway: 'inflammation', 'PCD', or 'PTM'\n"
           "  - Custom file path: '/path/to/your_pathway.csv'\n\n"
           "Example:\n"
           "  pathway_dict = scg.load_pathway_genes('inflammation')\n"
           "  pathway_dict = scg.load_pathway_genes('path/to/custom.csv')"
       )

   builtin_pathways = {
       'inflammation': 'inflammation.csv',
       'pcd': 'PCD.csv',
       'ptm': 'PTM.csv'
   }

   if pathway_file.lower() in builtin_pathways:
       pathway_name = pathway_file.lower()
       filename = builtin_pathways[pathway_name]

       try:
           import pkg_resources
           filepath = pkg_resources.resource_filename(
               'scGeneset', f'data/{filename}'
           )
       except ImportError:
           import importlib.resources as pkg_resources
           with pkg_resources.path('scGeneset.data', filename) as p:
               filepath = str(p)

       if not os.path.exists(filepath):
           raise FileNotFoundError(
               f"Error: Built-in pathway file '{filename}' not found!\n"
               f"Please check your scGeneset installation."
           )

       print(f"Loading built-in pathway: {pathway_name}")

   else:
       filepath = pathway_file

       if not os.path.exists(filepath):
           raise FileNotFoundError(
               f"Error: Custom pathway file not found: {filepath}\n"
               f"Please check the file path."
           )

       print(f"Loading custom pathway from: {os.path.basename(filepath)}")

   pathway_df = pd.read_csv(filepath)

   if pathway_df.shape[1] < 2:
       raise ValueError(
           f"Error: Invalid CSV format!\n"
           f"Expected at least 2 columns (pathway_name, gene_name)\n"
           f"Found {pathway_df.shape[1]} column(s)"
       )

   col_pathway = pathway_df.columns[0]
   col_gene = pathway_df.columns[1]

   pathway_dict = pathway_df.groupby(col_pathway)[col_gene].apply(set).to_dict()

   print(f"Loaded {len(pathway_dict)} pathways with {pathway_df.shape[0]} total genes")

   return pathway_dict
