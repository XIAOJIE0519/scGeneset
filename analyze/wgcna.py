"""
scGeneset/analyze/wgcna.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from ..utils.config import CMAP_POS
from ..utils.helper import simplify_name, get_residuals
from scipy.stats import spearmanr
import statsmodels.api as sm


def calc_wgcna(data, method='AUCell', group=None, cell_type=None, z=True,
              pathway=None, custom_genes=None):
   """
   Calculate WGCNA co-expression network between pathways

   Parameters:
       data: AnnData object
       method: Scoring method
       group: Group label (default: all cells)
       cell_type: Cell type (default: all cells)
       z: Whether to use normalized data
       pathway: Pathway list (if None, use all including Custom)
       custom_genes: Custom gene set

   Returns:
       fig: matplotlib Figure object
       tom_df: TOM matrix DataFrame
   """

   if group is None and cell_type is None:
       work_data = data.copy()
   elif group is not None and cell_type is not None:
       work_data = data[(data.obs['group'] == group) &
                        (data.obs['cell_type'] == cell_type)].copy()
   elif group is not None:
       work_data = data[data.obs['group'] == group].copy()
   else:
       work_data = data[data.obs['cell_type'] == cell_type].copy()

   if work_data.n_obs < 20:
       raise ValueError(f"Insufficient cell count ({work_data.n_obs})")

   suffix = '_Z' if z else ''
   if pathway is None:
       pathway_cols = [col for col in data.obs.columns
                       if col.startswith(f'{method}_') and col.endswith(suffix)]
   else:
       pathway_cols = [f'{method}_{p.replace("/", "_")}{suffix}' for p in pathway
                       if f'{method}_{p.replace("/", "_")}{suffix}' in data.obs.columns]

   if len(pathway_cols) < 2:
       raise ValueError("Insufficient number of pathways")

   pathway_scores = work_data.obs[pathway_cols].copy().dropna()

   if pathway_scores.empty:
       raise ValueError("Pathway score data is empty")

   scaler = StandardScaler()
   pathway_scores_scaled = pd.DataFrame(
       scaler.fit_transform(pathway_scores),
       index=pathway_scores.index,
       columns=pathway_scores.columns
   )

   print("Calculating adjusted correlation matrix...")
   corr_matrix = _calc_adjusted_corr_matrix(work_data, pathway_cols, method,
                                            pathway_scores_scaled)

   corr_matrix = corr_matrix.fillna(0)

   beta = 6
   adjacency = np.power(np.abs(corr_matrix.values), beta)
   adjacency_df = pd.DataFrame(adjacency, index=corr_matrix.index,
                               columns=corr_matrix.columns)

   print("Calculating TOM matrix...")
   tom = _calculate_tom(adjacency)
   tom_df = pd.DataFrame(tom, index=adjacency_df.index, columns=adjacency_df.columns)

   fig, ax = plt.subplots(figsize=(10, 9))

   sns.heatmap(tom_df, cmap=CMAP_POS, vmin=0, vmax=1, ax=ax, square=True,
               cbar_kws={'label': 'Topological Overlap', 'shrink': 0.8},
               xticklabels=[simplify_name(p.replace('_Z', ''), method) for p in pathway_cols],
               yticklabels=[simplify_name(p.replace('_Z', ''), method) for p in pathway_cols],
               linewidths=0.5, linecolor='white')

   z_label = " (Z-normalized)" if z else ""
   ax.set_title(f'WGCNA Co-expression Network (TOM Matrix)\n[{method}{z_label}]',
                fontsize=12, pad=10)
   plt.xticks(rotation=90, fontsize=9)
   plt.yticks(rotation=0, fontsize=9)
   plt.tight_layout()

   return fig, tom_df


def _calculate_tom(adj_matrix):
   """Calculate TOM matrix"""
   adj = adj_matrix.copy()
   np.fill_diagonal(adj, 0)
   k = adj.sum(axis=1)
   k_safe = np.where(k == 0, 1e-6, k)
   numerator = adj @ adj + adj
   denominator = np.minimum(k_safe[:, None], k_safe[None, :]) + 1 - adj
   denominator[denominator == 0] = 1e-6
   tom = numerator / denominator
   np.fill_diagonal(tom, 1)
   return tom


def _calc_adjusted_corr_matrix(data, pathway_cols, method, scores_scaled):
   """Calculate adjusted correlation matrix"""
   if 'pathway_dict' not in data.uns:
       print("Warning: pathway_dict not found in adata.uns, using empty dict")
       pathway_dict = {}
   else:
       print("Loaded pathway_dict from adata.uns:", len(data.uns['pathway_dict']), "pathways")
       pathway_dict = {k: set(v) for k, v in data.uns['pathway_dict'].items()}
   mat = data.raw.to_adata()
   gene_exp_array = mat.X.toarray() if hasattr(mat.X, 'toarray') else mat.X
   gene_names = mat.var_names.tolist()
   gene_to_idx = {g: i for i, g in enumerate(gene_names)}
   n_genes = mat.n_vars
   n_cells = data.n_obs

   n_pathways = len(pathway_cols)
   corr_matrix = pd.DataFrame(np.eye(n_pathways), index=pathway_cols,
                              columns=pathway_cols)

   for i in tqdm(range(n_pathways), desc="Calculating adjusted correlation"):
       for j in range(i + 1, n_pathways):
           pathway_A_col = pathway_cols[i]
           pathway_B_col = pathway_cols[j]

           raw_name_A = _get_raw_pathway_name(pathway_A_col, method, pathway_dict)
           raw_name_B = _get_raw_pathway_name(pathway_B_col, method, pathway_dict)

           if raw_name_A is None or raw_name_B is None:
               corr, _ = spearmanr(scores_scaled[pathway_A_col].dropna(),
                                   scores_scaled[pathway_B_col].dropna())
           else:
               genes_A = pathway_dict.get(raw_name_A, set())
               genes_B = pathway_dict.get(raw_name_B, set())
               genes_O = genes_A.intersection(genes_B)

               if len(genes_O) < 5:
                   corr, _ = spearmanr(scores_scaled[pathway_A_col].dropna(),
                                       scores_scaled[pathway_B_col].dropna())
               else:
                   indices_O = [gene_to_idx[g] for g in genes_O if g in gene_to_idx]
                   score_O = _calc_overlap_score_wgcna(indices_O, method,
                                                       gene_exp_array, n_cells,
                                                       n_genes, mat)

                   score_O_series = pd.Series(score_O, index=scores_scaled.index)
                   score_A = scores_scaled[pathway_A_col]
                   score_B = scores_scaled[pathway_B_col]

                   if score_O_series.std() > 0:
                       score_O_scaled = (score_O_series - score_O_series.mean()) / score_O_series.std()
                   else:
                       score_O_scaled = score_O_series

                   res_A = get_residuals(score_A, score_O_scaled)
                   res_B = get_residuals(score_B, score_O_scaled)

                   combined = pd.DataFrame({'A': res_A, 'B': res_B}).dropna()

                   if combined.shape[0] < 3:
                       corr = np.nan
                   else:
                       corr, _ = spearmanr(combined['A'], combined['B'])

           corr_matrix.iloc[i, j] = corr_matrix.iloc[j, i] = (
               corr if not np.isnan(corr) else 0.0
           )

   return corr_matrix


def _get_raw_pathway_name(pathway_col, method, pathway_dict):
   """Get raw pathway name"""
   prefix = f'{method}_'
   for key in pathway_dict.keys():
       safe_key = f"{prefix}{key.replace('/', '_')}"
       if pathway_col.startswith(safe_key):
           return key
   return None


def _calc_overlap_score_wgcna(indices, method, gene_exp_array, n_cells, n_genes, mat):
   """Calculate overlap gene score for WGCNA"""
   from ..core.scorer import (aucell_score_single_cell, ssgsea_score_single_cell,
                              calc_singscore_optimized, seurat_score_single)

   indices_np = np.array(indices)

   if method == 'AUCell':
       scores = [aucell_score_single_cell(i, gene_exp_array, indices_np)
                 for i in range(n_cells)]
   elif method == 'ssGSEA':
       scores = [ssgsea_score_single_cell(i, gene_exp_array, indices_np, n_genes)
                 for i in range(n_cells)]
   elif method == 'singscore':
       scores = [calc_singscore_optimized(i, gene_exp_array, indices_np, n_genes)
                 for i in range(n_cells)]
   else:
       scores = np.zeros(n_cells)

   return np.array(scores)
