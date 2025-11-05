"""
scGeneset/plot/corr.py (Correlation Analysis)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
import statsmodels.api as sm
from ..utils.config import CMAP_POS, GLOBAL_PALETTE
from ..utils.helper import simplify_name, calc_gene_score_subset, get_residuals
from ..core.loader import load_pathway_genes


def plot_corr_scatter(data, method='AUCell', group=None, cell_type=None,
                      z=True, custom_genes=None):
    """
    Plot scatter plot matrix for Custom vs 16 pathways (16 subplots)

    Parameters:
        data: AnnData object
        method: Scoring method
        group: Group (default None for all cells, e.g., 'HC', 'UC', 'CD')
        cell_type: Cell type (default None for all cells, e.g., 'Paneth cells')
        z: Whether to use normalized data
        custom_genes: Custom gene list

    Returns:
        fig: matplotlib Figure object
    """

    title_suffix = ""
    if group is None and cell_type is None:
        cell_adata = data.copy()
        title_suffix = "All Cells"
    elif group is not None and cell_type is not None:
        cell_adata = data[(data.obs['group'] == group) &
                          (data.obs['cell_type'] == cell_type)].copy()
        title_suffix = f"{cell_type} ({group})"
    elif group is not None:
        cell_adata = data[data.obs['group'] == group].copy()
        title_suffix = f"All Types ({group})"
    else:
        cell_adata = data[data.obs['cell_type'] == cell_type].copy()
        title_suffix = f"{cell_type} (All Groups)"

    if cell_adata.n_obs < 10:
        raise ValueError(f"Insufficient cells ({cell_adata.n_obs}) for {title_suffix}")

    gene_expr = calc_gene_score_subset(cell_adata, set(custom_genes), method=method)
    if len(gene_expr) == 0 or np.all(gene_expr == 0):
        raise ValueError(f"Custom gene set score is empty or all zeros for {title_suffix}")

    suffix = '_Z' if z else ''
    pathway_cols = [col for col in data.obs.columns
                    if col.startswith(f'{method}_') and col.endswith(suffix)
                    and 'Custom' not in col]

    if len(pathway_cols) == 0:
        raise ValueError("No valid pathways found matching criteria")

    n_scores = len(pathway_cols)
    n_cols = min(4, n_scores)
    n_rows = int(np.ceil(n_scores / n_cols))

    ridge_height = 0.6
    ridge_width = 0.6
    fig_width = 4 * n_cols
    fig_height = 4 * n_rows + ridge_height

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
    axes = axes.flatten() if n_scores > 1 else [axes]

    for i, score_col in enumerate(pathway_cols):
        ax = axes[i]
        pathway_scores = cell_adata.obs[score_col]
        color = GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]

        combined_df = pd.DataFrame({
            'gene_expr': gene_expr,
            'pathway_scores': pathway_scores
        }).dropna()

        if combined_df.empty or len(combined_df) < 3:
            corr, p_val = np.nan, np.nan
        else:
            corr, p_val = spearmanr(combined_df['gene_expr'],
                                    combined_df['pathway_scores'])

        sns.regplot(x=combined_df['gene_expr'], y=combined_df['pathway_scores'],
                    ax=ax, scatter_kws={'alpha': 0.1, 'color': 'grey', 's': 10},
                    line_kws={'color': color, 'linewidth': 2})

        base_color_rgb = mpatches.to_rgb(color) if isinstance(color, str) else color

        lighter_color = tuple(min(1.0, c + 0.3) for c in base_color_rgb[:3])
        top_ax = ax.inset_axes([0, 1, 1, ridge_height / 4], transform=ax.transAxes)
        sns.kdeplot(x=combined_df['gene_expr'], ax=top_ax, color=lighter_color,
                    fill=True, alpha=0.5)
        top_ax.set_xlim(ax.get_xlim())
        top_ax.axis('off')

        darker_color = tuple(max(0.0, c - 0.2) for c in base_color_rgb[:3])
        right_ax = ax.inset_axes([1, 0, ridge_width / 4, 1], transform=ax.transAxes)
        sns.kdeplot(y=combined_df['pathway_scores'], ax=right_ax, color=darker_color,
                    fill=True, alpha=0.5)
        right_ax.set_ylim(ax.get_ylim())
        right_ax.axis('off')

        cell_type = f"R = {corr:.2f}\np = {p_val:.2e}" if not np.isnan(corr) else "No data"
        ax.text(0.95, 0.05, cell_type, transform=ax.transAxes, ha='right',
                va='bottom', fontsize=9, bbox=dict(boxstyle='round,pad=0.3',
                                                   fc='white', alpha=0.7))

        ax.set_title(simplify_name(score_col.replace('_Z', ''), method), fontsize=11)
        ax.set_xlabel('Custom Expression', fontsize=10)
        ax.set_ylabel('Pathway Score', fontsize=10)
        sns.despine(ax=ax)
        ax.grid(False)

    for i in range(n_scores, len(axes)):
        axes[i].set_visible(False)

    z_label = " (Z-normalized)" if z else ""
    fig.suptitle(f'Custom vs Pathway Scores in {title_suffix}\n[{method}{z_label}]',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()

    return fig


def plot_corr_matrix(data, method='AUCell', group=None, cell_type=None,
                     z=True, custom_genes=None):
    """
    Plot correlation matrix for Custom vs pathways (half-matrix visualization with adjusted correlation)

    Parameters:
        data: AnnData object
        method: Scoring method
        group: Group (default None for all, e.g., 'HC', 'UC', 'CD')
        cell_type: Cell type (default None for all, e.g., 'Plasma cells')
        z: Whether to use normalized data
        custom_genes: Custom gene list

    Returns:
        fig: matplotlib Figure object
        corr_matrix: Correlation matrix DataFrame
    """

    if custom_genes is None or len(custom_genes) == 0:
        raise ValueError("custom_genes must be provided")

    title_suffix = ""
    if group is None and cell_type is None:
        cell_adata = data.copy()
        title_suffix = "All Cells"
    elif group is not None and cell_type is not None:
        cell_adata = data[(data.obs['group'] == group) &
                          (data.obs['cell_type'] == cell_type)].copy()
        title_suffix = f"{cell_type} ({group})"
    elif group is not None:
        cell_adata = data[data.obs['group'] == group].copy()
        title_suffix = f"All Types ({group})"
    else:
        cell_adata = data[data.obs['cell_type'] == cell_type].copy()
        title_suffix = f"{cell_type} (All Groups)"

    if cell_adata.n_obs < 10:
        raise ValueError(f"Insufficient cells ({cell_adata.n_obs}) for {title_suffix}")

    suffix = '_Z' if z else ''
    pathway_cols = [col for col in data.obs.columns
                    if col.startswith(f'{method}_') and col.endswith(suffix)
                    and 'Custom' not in col]

    n = len(pathway_cols)

    if n <= 5:
        fig_size = (8, 8)
    elif n <= 10:
        fig_size = (10, 10)
    elif n <= 20:
        fig_size = (14, 14)
    else:
        fig_size = (max(14, n * 0.7), max(14, n * 0.7))

    fig, ax = plt.subplots(figsize=fig_size)

    corr_matrix = pd.DataFrame(np.eye(n), index=pathway_cols, columns=pathway_cols)

    if 'pathway_dict' not in data.uns:
        raise ValueError(
            "Error: 'pathway_dict' not found in adata.uns!\n"
            "Please ensure you have run score_pathways() with pathway_dict specified."
        )

    pathway_dict = {k: set(v) for k, v in data.uns['pathway_dict'].items()}
    print(f"Loaded pathway_dict from adata.uns: {len(pathway_dict)} pathways")

    import anndata
    if data.raw is None:
        raise AttributeError("AnnData object must have a .raw attribute for gene expression data.")

    if not isinstance(data.raw, anndata.AnnData):
        if hasattr(data.raw, 'to_adata'):
            print("Warning: data.raw is not AnnData, converting...")
            data = data.copy()
            data.raw = data.raw.to_adata()
        else:
            raise TypeError(f"data.raw must be an AnnData object, but got {type(data.raw)}.")

    mat = data.raw.to_adata()
    gene_exp_array = mat.X.toarray() if hasattr(mat.X, 'toarray') else mat.X
    gene_names = mat.var_names.tolist()
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    n_genes = mat.n_vars

    cell_indices_global = np.where(data.obs.index.isin(cell_adata.obs.index))[0]
    sub_gene_exp_array = gene_exp_array[cell_indices_global, :]
    sub_n_cells = len(cell_indices_global)
    sub_adata_raw = mat[cell_indices_global, :].copy()

    scores_df = cell_adata.obs[pathway_cols]

    print(f"Calculating adjusted correlation matrix for {title_suffix}...")

    tasks = []
    for i in range(n):
        for j in range(i):
            tasks.append((i, j, pathway_cols[i], pathway_cols[j], pathway_dict,
                          method, gene_to_idx, sub_gene_exp_array, sub_n_cells,
                          n_genes, sub_adata_raw, scores_df, suffix))

    results = []
    for task_args in tasks:
        result = _calc_pathway_corr_task(*task_args)
        results.append(result)

    for (i, j, corr_val) in results:
        corr_matrix.iloc[i, j] = corr_val

    gene_expr = calc_gene_score_subset(cell_adata, set(custom_genes), method=method)
    gene_pathway_corr = {}

    if gene_expr is not None and not np.all(gene_expr == 0):
        for pathway_col in pathway_cols:
            combined_df = pd.DataFrame({
                'gene': gene_expr,
                'pathway': cell_adata.obs[pathway_col]
            }).dropna()

            if combined_df.shape[0] < 3:
                corr, pval = np.nan, 1.0
            else:
                corr, pval = spearmanr(combined_df['gene'], combined_df['pathway'])

            gene_pathway_corr[pathway_col] = {
                'corr': corr if not np.isnan(corr) else 0.0,
                'pval': pval
            }
    else:
        for pathway_col in pathway_cols:
            gene_pathway_corr[pathway_col] = {'corr': 0, 'pval': 1}

    fig, ax = plt.subplots(figsize=(10, 10))

    norm = plt.matplotlib.colors.Normalize(vmin=-1, vmax=1)
    cmap = plt.cm.RdBu_r

    for i in range(n):
        for j in range(i):
            corr_val = corr_matrix.iloc[i, j]
            abs_corr = abs(corr_val)
            size = 0.3 + 0.6 * abs_corr
            center_x = j + 0.5
            center_y = i + 0.5
            half_size = size / 2

            rect = plt.Rectangle(
                (center_x - half_size, center_y - half_size),
                size, size,
                facecolor=cmap(norm(corr_val)),
                edgecolor='black',
                linewidth=1.5,
                zorder=5
            )
            ax.add_patch(rect)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(0, n)
    ax.set_ylim(n, 0)
    ax.set_aspect('equal')
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_xticklabels([simplify_name(p.replace('_Z', ''), method) for p in pathway_cols],
                       rotation=90, fontsize=9)
    ax.set_yticklabels([simplify_name(p.replace('_Z', ''), method) for p in pathway_cols],
                       rotation=0, fontsize=9)
    ax.tick_params(axis='both', which='major', pad=0.1)
    z_label = " (Z-normalized)" if z else ""
    ax.set_title(f'Pathway Correlation in {title_suffix}\n[{method}{z_label}]',
                 fontsize=12, pad=10)
    ax.grid(False)

    custom_label_text = 'Custom'
    custom_pos_x = n - 5 + 0.5
    custom_pos_y = 4.5

    ax.text(custom_pos_x, custom_pos_y, custom_label_text,
            ha='center', va='center',
            fontweight='bold', color='black',
            bbox=dict(boxstyle='square,pad=0.3', fc='white', ec='none', alpha=0.9),
            zorder=20, transform=ax.transData)

    positive_corr_color = '#d62728'
    negative_corr_color = '#1f77b4'

    for i in range(n):
        corr_val = gene_pathway_corr[pathway_cols[i]]['corr']
        pval = gene_pathway_corr[pathway_cols[i]]['pval']

        color = positive_corr_color if corr_val > 0 else negative_corr_color if corr_val < 0 else 'grey'

        if pval < 0.05:
            abs_c = abs(corr_val)
            if abs_c <= 0.1:
                linewidth = 1
            elif abs_c <= 0.3:
                linewidth = 2
            elif abs_c <= 0.5:
                linewidth = 3
            else:
                linewidth = 4
        else:
            linewidth = 1.0

        final_color = color if pval < 0.05 else 'grey'

        start_pos = (custom_pos_x, custom_pos_y)
        end_pos = (i + 0.5, i + 0.5)

        rad = 0.1 if i < n / 2 else -0.1

        arc = mpatches.FancyArrowPatch(
            start_pos, end_pos,
            connectionstyle=f"arc3,rad={rad}",
            color=final_color,
            linewidth=linewidth,
            clip_on=False,
            zorder=10,
            arrowstyle='-|>', mutation_scale=10
        )
        ax.add_patch(arc)
        ax.plot(end_pos[0], end_pos[1], 'o', color='black', markersize=3, zorder=15)

    ax.set_xlim(-1.5, n + 1.5)
    ax.set_ylim(n + 1.5, -1.5)

    legend_handles = [
        Line2D([0], [0], color=positive_corr_color, lw=4, label='Positive Corr. (Custom vs Pathway, p<0.05)'),
        Line2D([0], [0], color=negative_corr_color, lw=4, label='Negative Corr. (Custom vs Pathway, p<0.05)'),
        Line2D([0], [0], color='grey', lw=2, label='Not Significant (p>=0.05)'),
        Line2D([0], [0], color='black', lw=1, label='|r| <= 0.1'),
        Line2D([0], [0], color='black', lw=2, label='0.1 < |r| <= 0.3'),
        Line2D([0], [0], color='black', lw=3, label='0.3 < |r| <= 0.5'),
        Line2D([0], [0], color='black', lw=4, label='|r| > 0.5')
    ]

    fig.legend(handles=legend_handles, loc='lower center',
               bbox_to_anchor=(0.5, -0.05), ncol=2, frameon=False, fontsize=8)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    fig.colorbar(sm, cax=cbar_ax, label='Correlation (r)')

    plt.subplots_adjust(left=0.12, right=0.9, bottom=0.18, top=0.92)

    return fig, corr_matrix


def _calc_pathway_corr_task(i, j, pathway_A_col, pathway_B_col, pathway_dict,
                            method, gene_to_idx, sub_gene_exp_array, sub_n_cells,
                            n_genes, sub_adata_raw, scores_df, suffix):
    """Task for calculating adjusted correlation between two pathways."""
    raw_name_A = None
    raw_name_B = None
    prefix = f'{method}_'

    for key in pathway_dict.keys():
        safe_key = f"{prefix}{key.replace('/', '_')}"
        if suffix:
            safe_key += suffix
        if safe_key == pathway_A_col:
            raw_name_A = key
        if safe_key == pathway_B_col:
            raw_name_B = key

    if raw_name_A is None or raw_name_B is None:
        corr, _ = spearmanr(scores_df[pathway_A_col].dropna(),
                            scores_df[pathway_B_col].dropna())
    else:
        genes_A = pathway_dict.get(raw_name_A, set())
        genes_B = pathway_dict.get(raw_name_B, set())
        genes_O = genes_A.intersection(genes_B)

        if len(genes_O) < 5:
            corr, _ = spearmanr(scores_df[pathway_A_col].dropna(),
                                scores_df[pathway_B_col].dropna())
        else:
            indices_O = [gene_to_idx[g] for g in genes_O if g in gene_to_idx]
            if not indices_O:
                corr, _ = spearmanr(scores_df[pathway_A_col].dropna(),
                                    scores_df[pathway_B_col].dropna())
            else:
                score_O = _calc_overlap_score(indices_O, method, sub_gene_exp_array,
                                              sub_n_cells, n_genes, sub_adata_raw)

                score_O_series = pd.Series(score_O, index=scores_df.index)
                score_A_series = scores_df[pathway_A_col]
                score_B_series = scores_df[pathway_B_col]

                if score_O_series.std() > 0:
                    score_O_scaled = (score_O_series - score_O_series.mean()) / score_O_series.std()
                else:
                    corr, _ = spearmanr(score_A_series.dropna(), score_B_series.dropna())
                    return (i, j, corr if not np.isnan(corr) else 0.0)

                res_A = get_residuals(score_A_series, score_O_scaled)
                res_B = get_residuals(score_B_series, score_O_scaled)

                combined_res = pd.DataFrame({'res_A': res_A, 'res_B': res_B}).dropna()

                if combined_res.shape[0] < 3:
                    corr = np.nan
                else:
                    corr, _ = spearmanr(combined_res['res_A'], combined_res['res_B'])

    return (i, j, corr if not np.isnan(corr) else 0.0)


def _calc_overlap_score(indices, method, gene_exp_array, n_cells, n_genes, adata_raw):
    """Calculate overlap gene score"""
    from ..core.scorer import (aucell_score_single_cell, ssgsea_score_single_cell,
                               calc_singscore_optimized, seurat_score_single)

    indices_np = np.array(indices)
    scores = np.zeros(n_cells)

    if method == 'AUCell':
        scores = [aucell_score_single_cell(cell_idx, gene_exp_array, indices_np)
                  for cell_idx in range(n_cells)]
    elif method == 'ssGSEA':
        scores = [ssgsea_score_single_cell(cell_idx, gene_exp_array, indices_np, n_genes)
                  for cell_idx in range(n_cells)]
    elif method == 'singscore':
        scores = [calc_singscore_optimized(cell_idx, gene_exp_array, indices_np, n_genes)
                  for cell_idx in range(n_cells)]
    elif method == 'Seurat':
        scores = [seurat_score_single(cell_idx, gene_exp_array, indices_np, n_genes)
                  for cell_idx in range(n_cells)]
    else:
        print(f"Warning: Scoring method '{method}' not handled. Returning zeros.")
        scores = np.zeros(n_cells)

    return np.array(scores)
