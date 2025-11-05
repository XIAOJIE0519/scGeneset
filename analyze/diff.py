"""
scGeneset/analyze/diff.py (Differential Analysis)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import matplotlib.patches as mpatches
from ..utils.helper import simplify_name


def analyze_differential(data, pathway, method='AUCell', group=None,
                        cell_type='Plasma cells', z=True):
   """
   Differential expression analysis (high score vs low score cells)

   Parameters:
       data: AnnData object
       pathway: Pathway name
       method: Scoring method
       group: Group label (default: all cells)
       cell_type: Cell type (default: 'Plasma cells')
       z: Whether to use normalized data

   Returns:
       results: Dictionary containing
           - 'fig_volcano': Volcano plot
           - 'deg_df': Differential gene DataFrame
   """

   suffix = '_Z' if z else ''
   pathway_col = f'{method}_{pathway.replace("/", "_")}{suffix}'

   if pathway_col not in data.obs.columns:
       raise ValueError(f"Pathway column not found: {pathway_col}")

   if group is not None:
       target_adata = data[(data.obs['cell_type'] == cell_type) &
                           (data.obs['group'] == group)].copy()
   else:
       target_adata = data[data.obs['cell_type'] == cell_type].copy()

   if target_adata.n_obs < 20:
       raise ValueError(f"Insufficient cell count ({target_adata.n_obs})")

   scores = target_adata.obs[pathway_col].dropna()
   score_quantile = scores.quantile([0.25, 0.75])
   low_cutoff = score_quantile.iloc[0]
   high_cutoff = score_quantile.iloc[1]

   target_adata.obs['score_group'] = 'Middle'
   target_adata.obs.loc[target_adata.obs[pathway_col] <= low_cutoff,
   'score_group'] = 'Low'
   target_adata.obs.loc[target_adata.obs[pathway_col] >= high_cutoff,
   'score_group'] = 'High'

   target_adata_filtered = target_adata[
       target_adata.obs['score_group'].isin(['Low', 'High'])
   ].copy()

   if len(target_adata_filtered[target_adata_filtered.obs['score_group'] == 'Low']) < 3:
       raise ValueError("Insufficient cells in Low group")
   if len(target_adata_filtered[target_adata_filtered.obs['score_group'] == 'High']) < 3:
       raise ValueError("Insufficient cells in High group")

   print("Running differential expression analysis...")
   sc.tl.rank_genes_groups(
       target_adata_filtered,
       groupby='score_group',
       groups=['High'],
       reference='Low',
       method='wilcoxon'
   )

   deg_df = sc.get.rank_genes_groups_df(target_adata_filtered, group='High')

   deg_df['-log10_pvals_adj'] = -np.log10(
       deg_df['pvals_adj'].replace(0, np.finfo(float).eps)
   )
   deg_df['logfoldchanges'] = deg_df['logfoldchanges'].fillna(0)

   deg_df['color'] = 'grey'
   lfc_threshold = 1
   pval_threshold = 0.05

   up_reg = (deg_df['logfoldchanges'] > lfc_threshold) & (deg_df['pvals_adj'] < pval_threshold)
   down_reg = (deg_df['logfoldchanges'] < -lfc_threshold) & (deg_df['pvals_adj'] < pval_threshold)
   deg_df.loc[up_reg, 'color'] = 'red'
   deg_df.loc[down_reg, 'color'] = 'blue'

   fig, ax = plt.subplots(figsize=(6, 5))

   sns.scatterplot(data=deg_df, x='logfoldchanges', y='-log10_pvals_adj',
                   hue='color', palette={'grey': 'grey', 'red': '#d62728',
                                         'blue': '#1f77b4'},
                   s=15, alpha=0.7, edgecolor='none', ax=ax)

   ax.axvline(x=lfc_threshold, linestyle='--', color='grey', linewidth=1)
   ax.axvline(x=-lfc_threshold, linestyle='--', color='grey', linewidth=1)
   ax.axhline(y=-np.log10(pval_threshold), linestyle='--', color='grey', linewidth=1)

   ax.set_xlim([-10, 10])
   ax.set_title(f'DEGs in High vs Low {simplify_name(pathway_col, method)} Score\n({cell_type})',
                fontsize=12)
   ax.set_xlabel('Log2 Fold Change', fontsize=10)
   ax.set_ylabel('-log10 (Adjusted P-value)', fontsize=10)

   legend_elements = [
       mpatches.Rectangle((0, 0), 1, 1, facecolor='#d62728', edgecolor='none',
                          label='Upregulated'),
       mpatches.Rectangle((0, 0), 1, 1, facecolor='#1f77b4', edgecolor='none',
                          label='Downregulated'),
       mpatches.Rectangle((0, 0), 1, 1, facecolor='grey', edgecolor='none',
                          label='Not Significant')
   ]
   ax.legend(handles=legend_elements, title='Regulation', frameon=False, fontsize=8)
   ax.grid(False)
   sns.despine(ax=ax)

   plt.tight_layout()

   results = {
       'fig_volcano': fig,
       'deg_df': deg_df
   }

   return results
