"""
scGeneset/analyze/pseudo.py - Pseudotime Analysis
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from sklearn.preprocessing import minmax_scale
import statsmodels.nonparametric.smoothers_lowess as lowess
from typing import List, Optional, Union

from ..utils.config import CMAP_POS, GLOBAL_PALETTE
from ..utils.helper import simplify_name


def analyze_pseudotime(data, method: str = 'AUCell', group: Optional[str] = None, z: bool = True,
                      pathway: Optional[List[str]] = None, custom_genes: Optional[dict] = None,
                      root_cell_type: str = 'Enterocytes',
                      root_group_for_pseudotime: Optional[str] = None):
   """
   Pseudotime analysis

   Parameters:
       data: AnnData object
       method: Scoring method
       group: The specific group to analyze. If None, analyze all cells.
              Note: If group is specified, root_group_for_pseudotime will be ignored,
              and root cells will be sought within this group.
       z: Whether to use Z-normalized pathway scores.
       pathway: List of pathway names to include. If None, use all available pathways.
       custom_genes: Custom gene sets
       root_cell_type: The cell type to use for identifying the pseudotime root.
                       Defaults to 'Enterocytes'.
       root_group_for_pseudotime: The group to use for identifying the pseudotime root.
                                  If group parameter is None, this parameter specifies
                                  which group to look for root_cell_type in.
                                  If None, and group is also None, the first available
                                  group from data.obs['group'] categories will be used.
                                  If group parameter is specified, this parameter is ignored.

   Returns:
       results: Dictionary containing:
           - 'fig_umap': Pseudotime UMAP figure
           - 'fig_heatmap': Heatmap figure
           - 'fig_trajectory': Trajectory plot figure
           - 'pseudotime': Pseudotime values as Series
   """

   if group is not None:
       work_data = data[data.obs['group'] == group].copy()
   else:
       work_data = data.copy()

   if work_data.n_obs < 50:
       raise ValueError(f"Insufficient cells ({work_data.n_obs}) in the selected subset for pseudotime analysis.")

   if 'X_umap' not in work_data.obsm:
       sc.tl.umap(work_data)

   sc.pp.neighbors(work_data, n_neighbors=15, use_rep='X_pca')
   sc.tl.diffmap(work_data)

   root_search_group = None
   if group is not None:
       root_search_group = group
   elif root_group_for_pseudotime is not None:
       root_search_group = root_group_for_pseudotime
   elif 'group' in data.obs and len(data.obs['group'].cat.categories) > 0:
       root_search_group = data.obs['group'].cat.categories[0]
       print(f"Info: No specific root group provided. Defaulting to '{root_search_group}' "
             f"for pseudotime root selection (using the first group found in data.obs['group']).")

   root_candidates_mask = (work_data.obs['cell_type'] == root_cell_type)
   if root_search_group is not None:
       root_candidates_mask = root_candidates_mask & (work_data.obs['group'] == root_search_group)

   root_candidates_indices = np.where(root_candidates_mask)[0]

   if len(root_candidates_indices) == 0:
       if root_search_group:
           print(f"Warning: No '{root_cell_type}' cells found in group '{root_search_group}' "
                 f"within the current subset. Using the first cell as pseudotime root.")
       else:
           print(f"Warning: No '{root_cell_type}' cells found within the current subset "
                 f"(or no group specified for root). Using the first cell as pseudotime root.")
       work_data.uns['iroot'] = 0
   else:
       work_data.uns['iroot'] = root_candidates_indices[0]

   sc.tl.dpt(work_data)
   work_data.obs['dpt_pseudotime_normalized'] = minmax_scale(
       work_data.obs['dpt_pseudotime']
   )

   fig_umap, ax = plt.subplots(figsize=(6, 5))
   fig_umap.patch.set_facecolor('white')
   ax.set_facecolor('white')

   z_label = " (Z-normalized)" if z else ""
   sc.pl.umap(work_data, color='dpt_pseudotime_normalized', cmap=CMAP_POS,
              size=20, title=f'Cell Trajectory by Pseudotime\n[{method}{z_label}]',
              ax=ax, show=False)

   plt.tight_layout()

   suffix = '_Z' if z else ''
   if pathway is None:
       pathway_cols = [col for col in data.obs.columns
                       if col.startswith(f'{method}_') and col.endswith(suffix)]
   else:
       pathway_cols = [f'{method}_{p.replace("/", "_")}{suffix}' for p in pathway
                       if f'{method}_{p.replace("/", "_")}{suffix}' in data.obs.columns]

   if len(pathway_cols) == 0:
       raise ValueError(f"No valid pathways found for method '{method}' {'(Z-normalized)' if z else ''}. "
                        "Please check pathway names or scoring method.")

   dpt_adata = work_data[~work_data.obs['dpt_pseudotime_normalized'].isna()].copy()

   if dpt_adata.n_obs == 0:
       raise ValueError("No valid pseudotime data available after filtering for NaNs.")

   pseudotime = dpt_adata.obs['dpt_pseudotime_normalized'].values

   bins = np.linspace(0, 1, 100)
   bin_indices = np.digitize(pseudotime, bins)

   heatmap_df = dpt_adata.obs[pathway_cols].copy()
   binned_scores = pd.DataFrame(np.nan, index=range(1, len(bins)),
                                columns=pathway_cols)

   for col in pathway_cols:
       for b in range(1, len(bins)):
           mask = bin_indices == b
           if np.any(mask):
               binned_scores.loc[b, col] = heatmap_df.loc[mask, col].mean()

   binned_scores_scaled = pd.DataFrame(
       minmax_scale(binned_scores.T, axis=1).T,
       index=binned_scores.index,
       columns=binned_scores.columns
   )

   fig_heatmap, ax = plt.subplots(figsize=(10, 4))

   sns.heatmap(binned_scores_scaled.T, cmap=CMAP_POS, ax=ax,
               xticklabels=False,
               yticklabels=[simplify_name(p.replace('_Z', ''), method) for p in pathway_cols],
               cbar_kws={'label': 'Normalized Score', 'shrink': 0.8,
                         'aspect': 30},
               linewidths=0)

   ax.set_xlabel('Pseudotime (0 to 1)', fontsize=10)
   ax.set_ylabel('Pathways', fontsize=10)
   z_label = " (Z-normalized)" if z else ""
   ax.set_title(f'Enrichment Trends along Pseudotime\n[{method}{z_label}]',
                fontsize=12, pad=10)

   x_ticks_positions = np.linspace(0, binned_scores_scaled.shape[0] - 1, 6)
   x_ticks_labels = ['{:.1f}'.format(x) for x in np.linspace(0, 1, 6)]
   ax.set_xticks(x_ticks_positions)
   ax.set_xticklabels(x_ticks_labels, rotation=0)

   plt.tight_layout()

   fig_trajectory, ax = plt.subplots(figsize=(14, 5))

   valid_indices = ~dpt_adata.obs['dpt_pseudotime_normalized'].isna()
   pseudotime_filtered = pseudotime[valid_indices.values]

   scores_df_filtered = dpt_adata.obs.loc[valid_indices, pathway_cols].copy()

   if scores_df_filtered.empty or scores_df_filtered.shape[0] < 2:
       print("Insufficient data for trajectory plot, skipping trajectory figure generation.")
       fig_trajectory = None
   else:
       sorted_idx = np.argsort(pseudotime_filtered)
       pseudotime_sorted = pseudotime_filtered[sorted_idx]

       for i, pathway_col in enumerate(pathway_cols):
           scores = scores_df_filtered[pathway_col].values[sorted_idx]

           valid_mask = ~np.isnan(scores)
           scores_valid = scores[valid_mask]
           pseudotime_valid = pseudotime_sorted[valid_mask]

           if len(scores_valid) > 1:
               smoothed = lowess.lowess(scores_valid, pseudotime_valid,
                                        frac=0.1, return_sorted=False)
               ax.plot(pseudotime_valid, smoothed,
                       label=simplify_name(pathway_col.replace('_Z', ''), method),
                       color=GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)],
                       linewidth=2.0, alpha=0.8)
           elif len(scores_valid) == 1:
               ax.plot(pseudotime_valid, scores_valid, 'o',
                       label=simplify_name(pathway_col.replace('_Z', ''), method),
                       color=GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)],
                       markersize=5)

       ax.set_xlabel('Pseudotime (0 to 1)', fontsize=10)
       ax.set_ylabel('Normalized Pathway Score', fontsize=10)
       z_label = " (Z-normalized)" if z else ""
       ax.set_title(f'Pathway Score Trajectories\n[{method}{z_label}]', fontsize=12)

       ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), ncol=1,
                 frameon=False, fontsize=8)
       sns.despine(ax=ax)
       ax.grid(False)
       plt.tight_layout(rect=[0, 0, 0.75, 1])

   results = {
       'fig_umap': fig_umap,
       'fig_heatmap': fig_heatmap,
       'fig_trajectory': fig_trajectory,
       'pseudotime': work_data.obs['dpt_pseudotime_normalized']
   }

   return results
