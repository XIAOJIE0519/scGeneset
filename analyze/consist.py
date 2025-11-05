"""
scGeneset/analyze/consist.py (Consistency Analysis)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from ..utils.config import CMAP_POS
from ..utils.helper import simplify_name


def calc_consistency(data, method='AUCell', group=None, z=True, custom_genes=None):
   """
   Calculate consistency between custom gene set and pathways across cell types

   Args:
       data: AnnData object
       method: Scoring method
       group: Group(s) to analyze (if None, calculate across all groups)
       z: Whether to use Z-normalized data
       custom_genes: Custom gene list

   Returns:
       fig: matplotlib Figure object
       consist_df: Consistency results DataFrame
   """

   if custom_genes is None or len(custom_genes) == 0:
       raise ValueError("custom_genes must be provided")

   import scanpy as sc
   gene_list = [g for g in custom_genes if g in data.raw.var_names]
   if not gene_list:
       raise ValueError("No genes from custom_genes found in data")

   sc.tl.score_genes(data, gene_list=gene_list, score_name='Custom_Signature',
                     use_raw=True)

   if group is not None:
       if isinstance(group, list):
           work_data = data[data.obs['group'].isin(group)].copy()
       else:
           work_data = data[data.obs['group'] == group].copy()
   else:
       work_data = data.copy()

   suffix = '_Z' if z else ''
   pathway_cols = [col for col in data.obs.columns
                   if col.startswith(f'{method}_') and col.endswith(suffix)
                   and 'Custom' not in col]

   consistency_results = []

   cell_types = sorted(work_data.obs['cell_type'].unique())
   cell_types_with_all = list(cell_types) + ['All Cells']

   for ct in cell_types_with_all:
       if ct == 'All Cells':
           adata_cell = work_data.copy()
       else:
           adata_cell = work_data[work_data.obs['cell_type'] == ct]

       if adata_cell.n_obs < 10:
           continue

       custom_scores = adata_cell.obs['Custom_Signature'].dropna()
       if custom_scores.empty:
           continue

       for pathway_col in pathway_cols:
           pathway_scores = adata_cell.obs[pathway_col].dropna()

           if pathway_scores.empty or len(custom_scores) < 3 or len(pathway_scores) < 3:
               continue

           common_index = custom_scores.index.intersection(pathway_scores.index)
           if len(common_index) < 3:
               continue

           corr, _ = spearmanr(custom_scores.loc[common_index],
                               pathway_scores.loc[common_index])

           if corr > 0.1:
               consistency_results.append({
                   'cell_type': ct,
                   'pathway': pathway_col,
                   'consistency': corr
               })

   consistency_df = pd.DataFrame(consistency_results)

   if consistency_df.empty:
       raise ValueError("No significant consistency found")

   fig, ax = plt.subplots(figsize=(10, 6))

   consistency_df['simplified_pathway'] = consistency_df['pathway'].apply(
       lambda p: simplify_name(p.replace('_Z', ''), method)
   )

   scatter = ax.scatter(
       x=consistency_df['cell_type'],
       y=consistency_df['simplified_pathway'],
       s=consistency_df['consistency'] * 500,
       c=consistency_df['consistency'],
       cmap=CMAP_POS,
       vmin=0.1,
       vmax=max(1.0, consistency_df['consistency'].max()),
       edgecolor='black',
       linewidth=0.5,
       alpha=0.8
   )

   z_label = " (Z-normalized)" if z else ""
   ax.set_title(f'Consistency: Custom Signature vs Pathways\n[{method}{z_label}]',
                fontsize=12, pad=10)
   ax.set_xlabel('Cell Subtype', fontsize=10)
   ax.set_ylabel('Pathway', fontsize=10)
   plt.xticks(rotation=45, ha='right', fontsize=8)
   plt.yticks(fontsize=8)
   sns.despine(ax=ax)
   ax.grid(False)

   cbar = plt.colorbar(scatter, ax=ax, pad=0.02, aspect=30)
   cbar.set_label('Correlation (r)', rotation=270, labelpad=10, fontsize=10)

   size_values = [0.2, 0.45, 0.8]
   size_labels = ["0.1 - 0.3", "0.3 - 0.6", "> 0.6"]
   size_legend_handles = [mpatches.Rectangle((0, 0), 1, 1, facecolor='none',
                                             edgecolor='none', label="< 0.1 (not shown)")]

   for i, val in enumerate(size_values):
       size_legend_handles.append(
           Line2D([0], [0], marker='o', color='w', markerfacecolor="gray",
                  markersize=np.sqrt(val * 500) / np.pi, label=size_labels[i],
                  linestyle='None')
       )

   fig.legend(handles=size_legend_handles, title='Correlation Size',
              loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=4,
              frameon=False, fontsize=8)

   plt.tight_layout(rect=[0, 0.1, 0.85, 1])

   return fig, consistency_df
