"""
scGeneset/plot/heat.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
from ..utils.config import CMAP_POS, GLOBAL_PALETTE
from ..utils.helper import simplify_name


def plot_heatmap(data, method='AUCell', group=None, cell_type=None, z=True,
                 pathway=None, custom_genes=None):
    """
    Plot pathway heatmap matrix.

    Args:
        data: AnnData object
        method: Scoring method
        group: Group list
        cell_type: Cell type list
        z: Whether to use Z-normalized data
        pathway: Pathway list
        custom_genes: Custom gene set

    Returns:
        fig: matplotlib Figure object
        data_matrix: Heatmap data matrix DataFrame
        p_matrix: P-value matrix DataFrame
    """

    suffix = '_Z' if z else ''
    if pathway is None:
        pathway_cols = [col for col in data.obs.columns
                        if col.startswith(f'{method}_') and col.endswith(suffix)]
    else:
        pathway_cols = [f'{method}_{p.replace("/", "_")}{suffix}' for p in pathway
                        if f'{method}_{p.replace("/", "_")}{suffix}' in data.obs.columns]

    if len(pathway_cols) == 0:
        raise ValueError(f"No valid pathway score columns found (method={method}, z={z})")

    if cell_type is None:
        cell_types = sorted(data.obs['cell_type'].unique())
        cell_types = list(cell_types) + ['All Cells']
    else:
        cell_types = cell_type if isinstance(cell_type, list) else [cell_type]

    if group is None:
        diagnoses = sorted(data.obs['group'].unique())
    else:
        diagnoses = group if isinstance(group, list) else [group]

    n_pathways = len(pathway_cols)
    n_groups = len(diagnoses)
    n_types = len(cell_types)

    data_matrix = np.zeros((n_types, n_pathways * n_groups))
    p_matrix = np.ones((n_types, n_pathways))

    for i, ct in enumerate(cell_types):
        if ct == 'All Cells':
            plot_df = data.obs.copy()
        else:
            plot_df = data.obs[data.obs['cell_type'] == ct].copy()

        if plot_df.empty:
            continue

        for j, pathway_col in enumerate(pathway_cols):
            group_scores = []
            for k, diag in enumerate(diagnoses):
                scores = plot_df[plot_df['group'] == diag][pathway_col].dropna()
                mean_score = scores.mean() if len(scores) > 0 else np.nan
                data_matrix[i, j * n_groups + k] = mean_score
                group_scores.append(scores)

            valid_groups = [g for g in group_scores if len(g) > 3]
            if len(valid_groups) >= 2:
                _, p_val = kruskal(*valid_groups)
                p_matrix[i, j] = p_val

    df_heatmap = pd.DataFrame(
        data_matrix,
        index=cell_types,
        columns=[f'{pathway_cols[i]}_{g}' for i in range(n_pathways) for g in diagnoses]
    )

    num_cols = n_pathways * n_groups
    num_rows = n_types
    cell_size = 0.5
    fig_width = num_cols * cell_size + 2
    fig_height = num_rows * cell_size + 2

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    mask = df_heatmap.isna()

    sns.heatmap(df_heatmap, cmap=CMAP_POS, ax=ax, annot=False,
                linewidths=0, square=True, mask=mask,
                cbar_kws={'label': 'Mean Pathway Score', 'shrink': 0.8,
                          'aspect': 30, 'pad': 0.02})

    for i in range(p_matrix.shape[0]):
        for j in range(p_matrix.shape[1]):
            if p_matrix[i, j] > 0.05:
                for k in range(n_groups):
                    x_pos = j * n_groups + k + 0.5
                    y_pos = i + 0.5
                    ax.text(x_pos, y_pos, 'X', ha='center', va='center',
                            color='gray', fontsize=10, fontweight='bold')

    for i in range(num_rows + 1):
        ax.axhline(i, color='black', linewidth=0.8)
    for j in range(0, num_cols + 1, n_groups):
        ax.axvline(j, color='black', linewidth=0.8)

    pathway_centers = [n_groups * idx + n_groups / 2 for idx in range(n_pathways)]
    ax.set_xticks(pathway_centers)
    ax.set_xticklabels([simplify_name(p.replace('_Z', ''), method) for p in pathway_cols],
                       rotation=90, ha='center', fontsize=8)
    ax.set_yticks([i + 0.5 for i in range(num_rows)])
    ax.set_yticklabels(cell_types, rotation=0, fontsize=8)
    ax.set_xlabel('')
    ax.set_ylabel('Cell Types', fontsize=10)
    z_label = " (Z-normalized)" if z else ""
    ax.set_title(f'Pathway Scores by Cell Type and Group\n[{method}{z_label}]',
                 fontsize=12, pad=10)

    group_categories = sorted(data.obs['group'].unique())
    current_group_colors = {
        g: GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
        for i, g in enumerate(group_categories)
    }

    for j in range(num_cols):
        g_idx = j % n_groups
        g = diagnoses[g_idx]
        color = current_group_colors.get(g, '#CCCCCC')
        ax.add_patch(plt.Rectangle((j, -0.4), 1, 0.4, color=color,
                                   clip_on=False, transform=ax.transData, zorder=5))

    plt.tight_layout()

    return fig, df_heatmap, pd.DataFrame(p_matrix, index=cell_types, columns=pathway_cols)
