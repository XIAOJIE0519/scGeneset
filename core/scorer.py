"""
scGeneset/core/scorer.py
"""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata
from tqdm import tqdm
from sklearn.preprocessing import minmax_scale


def aucell_score_single_cell(cell_idx, gene_exp_array, gene_indices):
   """AUCell method for single-cell scoring"""
   expression = gene_exp_array[cell_idx]
   if len(gene_indices) == 0:
       return 0.0
   sorted_indices = np.argsort(expression)[::-1]
   ranks = np.zeros(len(expression))
   ranks[sorted_indices] = np.arange(len(expression))
   gene_ranks = ranks[gene_indices]
   max_rank = len(expression)
   auc = 1 - (np.mean(gene_ranks) / max_rank)
   return auc


def ssgsea_score_single_cell(cell_idx, gene_exp_array, gene_indices, n_genes):
   """ssGSEA method for single-cell scoring"""
   expression = gene_exp_array[cell_idx]
   if len(gene_indices) == 0:
       return 0.0
   sorted_indices = np.argsort(expression)[::-1]
   rank_vector = np.zeros(n_genes)
   rank_vector[gene_indices] = 1
   sorted_rank = rank_vector[sorted_indices]
   cumsum = np.cumsum(sorted_rank)
   es = np.sum(cumsum) / (len(gene_indices) * n_genes)
   return es


def calc_singscore_optimized(cell_idx, gene_exp_array, gene_set_indices, n_genes):
   """singscore method for single-cell scoring"""
   expr = gene_exp_array[cell_idx]
   ranks = rankdata(expr) / n_genes
   if len(gene_set_indices) > 0:
       score = np.mean(ranks[gene_set_indices])
   else:
       score = 0.0
   return score


def seurat_score_single(pathway, genes, adata_raw, ctrl_size=100):
   """AddModuleScore (Seurat-style)"""
   genes_in_data = [g for g in genes if g in adata_raw.var_names]
   if len(genes_in_data) < 2:
       return pathway, np.zeros(adata_raw.n_obs)
   temp_adata = adata_raw.copy()
   sc.tl.score_genes(
       temp_adata,
       genes_in_data,
       ctrl_size=min(ctrl_size, len(genes_in_data)),
       score_name='temp_score',
       use_raw=False
   )
   scores = temp_adata.obs['temp_score'].values
   return pathway, scores


def score_pathways(adata, pathway_dict=None, custom_genes=None,
                  methods=['AUCell', 'AddModule', 'ssGSEA', 'singscore'],
                  normalize=True, min_genes=2):
   """
   Calculate pathway activity scores

   Parameters:
       adata: AnnData object
       pathway_dict: {pathway_name: gene_set} dictionary, auto-loaded if None
       custom_genes: Custom gene list
       methods: List of scoring methods
       normalize: Whether to normalize (z=True)
       min_genes: Minimum gene count threshold

   Returns:
       adata: AnnData object with added scoring results
       method_stats: Statistics dictionary
   """
   print("=" * 60)
   print("Starting pathway scoring")
   print("=" * 60)

   if pathway_dict is None:
       from .loader import load_pathway_genes
       raise ValueError(
           "Error: 'pathway_dict' must be provided!\n"
           "Please load pathways first:\n"
           "  pathway_dict = scg.load_pathway_genes('inflammation')\n"
           "Then pass it to score_pathways:\n"
           "  adata, stats = scg.score_pathways(adata, pathway_dict=pathway_dict)"
       )

   if custom_genes is not None and len(custom_genes) > 0:
       pathway_dict['Custom'] = set(custom_genes)
       print(f"Added custom gene set 'Custom': {len(custom_genes)} genes")

   pathway_dict_serializable = {}
   for pathway_name, genes in pathway_dict.items():
       safe_name = pathway_name.replace('/', '_')
       pathway_dict_serializable[safe_name] = list(genes)
   adata.uns['pathway_dict'] = pathway_dict_serializable
   print("Saved pathway_dict to adata.uns['pathway_dict']")

   mat = adata.raw.to_adata()
   gene_exp_array = mat.X.toarray() if hasattr(mat.X, 'toarray') else mat.X
   gene_names = mat.var_names.tolist()
   gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}

   pathway_gene_indices = {}
   for pathway, genes in pathway_dict.items():
       indices = [gene_to_idx[g] for g in genes if g in gene_to_idx]
       if len(indices) >= min_genes:
           pathway_gene_indices[pathway] = np.array(indices)
       else:
           print(f"Skipping '{pathway}': only {len(indices)} genes found (requires >={min_genes})")

   print(f"Scoring {len(pathway_gene_indices)} pathways")

   method_stats = {}

   for method in methods:
       print(f"\n{'─' * 60}")
       print(f"Method: {method}")
       print(f"{'─' * 60}")

       if method == 'AUCell':
           results = [
               _compute_aucell(pathway, indices, gene_exp_array, adata.n_obs)
               for pathway, indices in tqdm(pathway_gene_indices.items(), desc=method)
           ]
       elif method == 'AddModule':
           results = [
               seurat_score_single(pathway, list(genes), mat)
               for pathway, genes in tqdm(pathway_dict.items(), desc=method)
               if pathway in pathway_gene_indices
           ]
       elif method == 'ssGSEA':
           results = [
               _compute_ssgsea(pathway, indices, gene_exp_array,
                                        adata.n_obs, mat.n_vars)
               for pathway, indices in tqdm(pathway_gene_indices.items(), desc=method)
           ]
       elif method == 'singscore':
           results = [
               _compute_singscore(pathway, indices, gene_exp_array,
                                           adata.n_obs, mat.n_vars)
               for pathway, indices in tqdm(pathway_gene_indices.items(), desc=method)
           ]
       else:
           print(f"Unknown method: {method}")
           continue

       scores_df = pd.DataFrame(index=adata.obs.index)
       added_cols = []

       for pathway, scores in results:
           pathway_safe = pathway.replace("/", "_")
           col_name = f'{method}_{pathway_safe}'
           scores_df[col_name] = scores
           added_cols.append(col_name)

           if normalize:
               col_name_z = f'{method}_{pathway_safe}_Z'
               scores_df[col_name_z] = minmax_scale(scores)
               added_cols.append(col_name_z)

       adata.obs = pd.concat([adata.obs, scores_df], axis=1)

       method_stats[method] = {
           'n_pathways': len(results),
           'n_columns': len(added_cols),
           'columns': added_cols
       }

       print(f"Completed {method}: added {len(scores_df.columns)} columns")

   print("\n" + "=" * 60)
   print(f"Scoring completed! Calculated {len(methods)} methods")
   print("=" * 60)

   return adata, method_stats


def _compute_aucell(pathway, gene_indices, gene_exp_array, n_cells):
   """Compute AUCell pathway scores"""
   scores = [
       aucell_score_single_cell(i, gene_exp_array, gene_indices)
       for i in range(n_cells)
   ]
   return pathway, scores


def _compute_ssgsea(pathway, gene_indices, gene_exp_array, n_cells, n_genes):
   """Compute ssGSEA pathway scores"""
   scores = [
       ssgsea_score_single_cell(i, gene_exp_array, gene_indices, n_genes)
       for i in range(n_cells)
   ]
   return pathway, scores


def _compute_singscore(pathway, gene_indices, gene_exp_array, n_cells, n_genes):
   """Compute singscore pathway scores"""
   scores = [
       calc_singscore_optimized(i, gene_exp_array, gene_indices, n_genes)
       for i in range(n_cells)
   ]
   return pathway, scores
