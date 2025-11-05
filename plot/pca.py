"""
scGeneset/plot/pca.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.spatial import ConvexHull
from scipy.interpolate import splprep, splev
from sklearn.preprocessing import minmax_scale
import matplotlib.patches as mpatches
from ..utils.config import CMAP_POS, GLOBAL_PALETTE
from ..utils.helper import simplify_name


def plot_pca(data, method='AUCell', z=True, pathway=None):
    """
    Plot PCA scatter plot and feature importance plot

    Parameters:
        data: AnnData object
        method: Scoring method
        z: Whether to use normalized data
        pathway: List of pathways (if None, use all)

    Returns:
        fig_scatter: PCA scatter plot
        fig_importance: Feature importance plot
        loadings_df: Loadings DataFrame
    """

    suffix = '_Z' if z else ''
    if pathway is None:
        pathway_cols = [col for col in data.obs.columns
                        if col.startswith(f'{method}_') and col.endswith(suffix)]
    else:
        pathway_cols = [f'{method}_{p.replace("/", "_")}{suffix}' for p in pathway
                        if f'{method}_{p.replace("/", "_")}{suffix}' in data.obs.columns]

    if len(pathway_cols) < 2:
        raise ValueError(f"Insufficient number of pathways ({len(pathway_cols)}), at least 2 required")

    pca_data = data.obs[pathway_cols].dropna(axis=1, how='all')

    if pca_data.shape[1] < 2:
        raise ValueError("Insufficient valid pathways")

    import anndata
    adata_pca = anndata.AnnData(pca_data)
    adata_pca.obs['group'] = data.obs['group'].values
    adata_pca.obs['group'] = adata_pca.obs['group'].astype('category')

    sc.tl.pca(adata_pca)

    group_categories = sorted(adata_pca.obs['group'].cat.categories)
    current_group_colors = {
        g: GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
        for i, g in enumerate(group_categories)
    }

    fig_scatter, ax = plt.subplots(figsize=(6, 5))

    for g in adata_pca.obs['group'].cat.categories:
        mask = adata_pca.obs['group'] == g
        coords = adata_pca.obsm['X_pca'][mask, :2]
        ax.scatter(coords[:, 0], coords[:, 1],
                   color=current_group_colors.get(g, '#CCCCCC'),
                   s=30, alpha=0.2, label=g, edgecolors='none')

    for g in adata_pca.obs['group'].cat.categories:
        points = adata_pca.obsm['X_pca'][adata_pca.obs['group'] == g, :2]
        if len(points) > 2:
            try:
                hull = ConvexHull(points)
                x_hull = points[hull.vertices, 0]
                y_hull = points[hull.vertices, 1]
                tck, u = splprep([x_hull, y_hull], s=0, per=True)
                u_new = np.linspace(u.min(), u.max(), 1000)
                x_new, y_new = splev(u_new, tck, der=0)
                ax.plot(x_new, y_new, color=current_group_colors.get(g, '#CCCCCC'),
                        lw=1.5, alpha=0.8)
            except:
                pass

    variance_ratio = adata_pca.uns['pca']['variance_ratio']
    ax.set_xlabel(f'PC1 ({variance_ratio[0] * 100:.2f}%)', fontsize=10)
    ax.set_ylabel(f'PC2 ({variance_ratio[1] * 100:.2f}%)', fontsize=10)
    z_label = " (Z-normalized)" if z else ""
    ax.set_title(f'PCA of Pathway Scores\n[{method}{z_label}]', fontsize=12, pad=10)

    group_counts = adata_pca.obs['group'].value_counts()

    legend_elements = [
        mpatches.Rectangle((0, 0), 1, 1,
                           facecolor=current_group_colors.get(g, '#CCCCCC'),
                           edgecolor='none',
                           label=f'{g} (n={group_counts.get(g, 0)})')
        for g in group_categories
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=8)

    sns.despine(ax=ax)
    ax.grid(False)
    plt.tight_layout()

    loadings = adata_pca.varm['PCs'][:, :2]
    loadings_df = pd.DataFrame(loadings, index=adata_pca.var_names,
                               columns=['PC1', 'PC2'])
    loadings_df.index = [simplify_name(idx.replace('_Z', ''), method) for idx in loadings_df.index]
    loadings_df['Total_Contribution'] = loadings_df[['PC1', 'PC2']].abs().sum(axis=1)
    loadings_df = loadings_df.sort_values('Total_Contribution', ascending=False)

    fig_importance, ax = plt.subplots(figsize=(8, 5))

    contrib_norm = minmax_scale(loadings_df['Total_Contribution'])
    colors = CMAP_POS(contrib_norm)

    sns.barplot(x='Total_Contribution', y=loadings_df.index, data=loadings_df,
                palette=list(colors), ax=ax)

    z_label = " (Z-normalized)" if z else ""
    ax.set_title(f'PCA Feature Importance\n[{method}{z_label}]', fontsize=12)
    ax.set_xlabel('Total Contribution (PC1 + PC2)', fontsize=10)
    ax.set_ylabel('Pathway', fontsize=10)
    sns.despine(ax=ax)
    ax.grid(False)

    import matplotlib.colors as mcolors
    cbar = plt.colorbar(
        plt.cm.ScalarMappable(
            cmap=CMAP_POS,
            norm=mcolors.Normalize(vmin=loadings_df['Total_Contribution'].min(),
                                   vmax=loadings_df['Total_Contribution'].max())
        ),
        ax=ax, pad=0.02, aspect=30
    )
    cbar.set_label('Contribution', rotation=270, labelpad=10)

    plt.tight_layout()

    return fig_scatter, fig_importance, loadings_df
