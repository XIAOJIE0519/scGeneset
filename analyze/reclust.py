"""
scGeneset/analyze/reclust.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from ..utils.config import GLOBAL_PALETTE, CMAP_POS
from ..utils.helper import simplify_name


def recluster_cells(data, pathway, method='AUCell', group=None, cell_type='Plasma cells',
                   z=True, cluster='high'):
   """
   Recluster cells with high or low scores for a specific pathway.

   Args:
       data: AnnData object
       pathway: Pathway name
       method: Scoring method
       group: Group filter (default: all cells)
       cell_type: Cell type (default: 'Plasma cells')
       z: Whether to use Z-normalized data
       cluster: 'high' or 'low', select high or low scoring cells

   Returns:
       results: Dictionary containing:
           - 'fig_umap': UMAP plot
           - 'fig_bubble': Marker gene bubble plot
           - 'adata_sub': Reclustered AnnData object
           - 'hvgs': Highly variable genes DataFrame
           - 'markers': Marker genes DataFrame
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

   if target_adata.n_obs < 30:
       raise ValueError(f"Insufficient number of cells ({target_adata.n_obs})")

   scores = target_adata.obs[pathway_col].dropna()

   if cluster == 'high':
       cutoff = scores.quantile(0.5)
       selected_cells = target_adata[target_adata.obs[pathway_col] >= cutoff].copy()
       title_suffix = "High"
   elif cluster == 'low':
       cutoff = scores.quantile(0.5)
       selected_cells = target_adata[target_adata.obs[pathway_col] < cutoff].copy()
       title_suffix = "Low"
   else:
       raise ValueError("cluster must be 'high' or 'low'")

   print(f"Selected {selected_cells.n_obs} {title_suffix}-scoring cells")

   if selected_cells.n_obs < 30:
       raise ValueError(f"Insufficient number of {title_suffix}-scoring cells")

   sc.pp.highly_variable_genes(selected_cells, n_top_genes=2000)
   sc.tl.pca(selected_cells, n_comps=30)
   sc.pp.neighbors(selected_cells, n_neighbors=15, n_pcs=30)
   sc.tl.leiden(selected_cells, resolution=0.5, key_added='subcluster')
   sc.tl.umap(selected_cells)

   cluster_counts = selected_cells.obs['subcluster'].value_counts()
   top5_clusters = cluster_counts.head(5).index.tolist()
   n_clusters_display = len(top5_clusters)

   print(f"Total {len(cluster_counts)} subclusters, displaying top {n_clusters_display} by cell count")
   for cluster_id in top5_clusters:
       print(f"  - Subcluster {cluster_id}: {cluster_counts[cluster_id]} cells")

   subcluster_hvgs = {}
   for cluster_id in selected_cells.obs['subcluster'].unique():
       cluster_cells = selected_cells[selected_cells.obs['subcluster'] == cluster_id]
       if cluster_cells.n_obs >= 10:
           sc.pp.highly_variable_genes(cluster_cells, n_top_genes=5, subset=False)
           hvg_genes = cluster_cells.var_names[cluster_cells.var['highly_variable']][:5].tolist()
           subcluster_hvgs[f'Subcluster_{cluster_id}'] = hvg_genes

   hvg_df = pd.DataFrame.from_dict(subcluster_hvgs, orient='index').T

   fig_umap, ax = plt.subplots(figsize=(8, 8))

   n_clusters = len(selected_cells.obs['subcluster'].unique())
   subcluster_colors = {str(i): GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
                        for i in range(n_clusters)}

   for cluster_id in selected_cells.obs['subcluster'].unique():
       mask = selected_cells.obs['subcluster'] == cluster_id
       cluster_coords = selected_cells.obsm['X_umap'][mask]
       cluster_count = mask.sum()
       ax.scatter(cluster_coords[:, 0], cluster_coords[:, 1],
                  color=subcluster_colors[str(cluster_id)], s=20, alpha=0.6,
                  label=f'Subcluster {cluster_id} ({cluster_count})',
                  edgecolors='none')

   sns.kdeplot(x=selected_cells.obsm['X_umap'][:, 0],
               y=selected_cells.obsm['X_umap'][:, 1],
               ax=ax, color='gray', levels=5, linewidths=1.5,
               alpha=0.7, linestyles='dashed')

   ax.set_xlabel('UMAP_1', fontsize=10)
   ax.set_ylabel('UMAP_2', fontsize=10)
   z_label = " (Z-normalized)" if z else ""
   ax.set_title(
       f'{cell_type} Subclusters ({title_suffix} {simplify_name(pathway_col.replace("_Z", ""), method)})\n[{method}{z_label}]',
       fontsize=12, pad=10)
   ax.legend(loc='best', frameon=False, fontsize=8)
   ax.set_aspect('equal')
   plt.tight_layout()

   print("Finding marker genes...")
   sc.tl.rank_genes_groups(selected_cells, groupby='subcluster', method='wilcoxon')

   top_genes_dict = {}
   for cluster_id in top5_clusters:
       cluster_genes = sc.get.rank_genes_groups_df(selected_cells, group=cluster_id)
       sig_genes = cluster_genes[
           (cluster_genes['pvals_adj'] < 0.05) &
           (cluster_genes['logfoldchanges'] > 0.5)
       ].head(5)
       top_genes_dict[f'Subcluster_{cluster_id}'] = sig_genes['names'].tolist()

   all_marker_genes = []
   for cluster_name, genes in top_genes_dict.items():
       for gene in genes:
           all_marker_genes.append((gene, cluster_name))

   if len(all_marker_genes) > 0:
       bubble_data = []

       selected_cells_top5 = selected_cells[selected_cells.obs['subcluster'].isin(top5_clusters)].copy()

       for gene, source_cluster in all_marker_genes:
           if gene not in selected_cells_top5.raw.var_names:
               continue

           gene_expr = selected_cells_top5.raw[:, gene].X.toarray().flatten()

           gene_label = f"{gene} (C{source_cluster.split('_')[1]})"

           cluster_num = source_cluster.split('_')[1]
           gene_color = subcluster_colors[cluster_num]

           for cluster_id in top5_clusters:
               cluster_mask = selected_cells_top5.obs['subcluster'] == cluster_id
               cluster_expr = gene_expr[cluster_mask]

               if len(cluster_expr) > 0:
                   mean_expr = np.mean(cluster_expr)
                   bubble_data.append({
                       'gene_label': gene_label,
                       'gene': gene,
                       'source_cluster': source_cluster,
                       'gene_color': gene_color,
                       'subcluster': f'Subcluster_{cluster_id}',
                       'expression': mean_expr
                   })

       bubble_df = pd.DataFrame(bubble_data)

       if not bubble_df.empty:
           from sklearn.preprocessing import minmax_scale
           bubble_df['expression_scaled'] = minmax_scale(bubble_df['expression'])

           unique_gene_labels = list(dict.fromkeys(bubble_df['gene_label']))

           gene_label_colors = {}
           for gene_label in unique_gene_labels:
               gene_label_colors[gene_label] = bubble_df[bubble_df['gene_label'] == gene_label]['gene_color'].iloc[0]

           fig_bubble, ax = plt.subplots(figsize=(max(8, n_clusters_display * 1.5),
                                                   max(6, len(unique_gene_labels) * 0.4)))

           scatter = ax.scatter(
               x=bubble_df['subcluster'],
               y=bubble_df['gene_label'],
               s=bubble_df['expression_scaled'] * 500 + 50,
               c=bubble_df['expression_scaled'],
               cmap=CMAP_POS,
               vmin=0,
               vmax=1,
               edgecolor='black',
               linewidth=0.5,
               alpha=0.8
           )

           ax.set_yticks(range(len(unique_gene_labels)))
           ax.set_yticklabels(unique_gene_labels, fontsize=9)

           for i, gene_label in enumerate(unique_gene_labels):
               ax.get_yticklabels()[i].set_color(gene_label_colors[gene_label])
               ax.get_yticklabels()[i].set_fontweight('bold')

           ax.set_xlabel('Subcluster', fontsize=10)
           ax.set_ylabel('Marker Genes', fontsize=10)
           z_label = " (Z-normalized)" if z else ""
           ax.set_title(f'Top Marker Genes per Subcluster\n[{method}{z_label}]',
                        fontsize=12, pad=10, fontweight='bold')

           plt.xticks(rotation=45, ha='right', fontsize=9)
           sns.despine(ax=ax)
           ax.grid(False)

           cbar = plt.colorbar(scatter, ax=ax, pad=0.02, aspect=30)
           cbar.set_label('Normalized Expression', rotation=270, labelpad=15, fontsize=10)

           size_values = [0.2, 0.5, 0.8]
           size_labels = ["Low", "Medium", "High"]
           size_legend_handles = []

           for i, val in enumerate(size_values):
               size_legend_handles.append(
                   Line2D([0], [0], marker='o', color='w',
                         markerfacecolor="gray",
                         markersize=np.sqrt(val * 500 + 50) / 2,
                         label=size_labels[i],
                         linestyle='None',
                         markeredgecolor='black',
                         markeredgewidth=0.5)
               )

           ax.legend(handles=size_legend_handles,
                    title='Expression Level',
                    loc='upper left',
                    bbox_to_anchor=(1.15, 1),
                    frameon=False,
                    fontsize=8)

           plt.tight_layout()

           marker_df = pd.DataFrame.from_dict(top_genes_dict, orient='index').T
       else:
           print("Warning: No valid expression data found")
           fig_bubble = None
           marker_df = pd.DataFrame()
   else:
       print("Warning: No significant marker genes found")
       fig_bubble = None
       marker_df = pd.DataFrame()

   results = {
       'fig_umap': fig_umap,
       'fig_bubble': fig_bubble,
       'adata_sub': selected_cells,
       'hvgs': hvg_df,
       'markers': marker_df
   }

   return results
