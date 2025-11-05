"""
scGeneset/plot/violin.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import matplotlib.patches as mpatches
from ..utils.config import GLOBAL_PALETTE
from ..utils.helper import simplify_name, format_pval


def plot_violin(data, pathway, method='AUCell', group=None, cell_type=None,
                z=True, custom_genes=None):
    """
    Plot pathway violin plot

    Parameters:
        data: AnnData object
        pathway: Pathway name or 'Custom'
        method: Scoring method
        group: If specified, compare different types within this group.
               Can be a single group string or a list of groups.
        cell_type: If specified, compare different groups within this cell type.
                   Can be a single cell type string or a list of cell types.
        z: Whether to use Z-normalized data
        custom_genes: Custom gene set

    Returns:
        fig: matplotlib Figure object
    """

    suffix = '_Z' if z else ''
    if pathway == 'Custom' or (custom_genes is not None):
        pathway_col = f'{method}_Custom{suffix}'
    else:
        pathway_col = f'{method}_{pathway.replace("/", "_")}{suffix}'

    if pathway_col not in data.obs.columns:
        raise ValueError(f"Pathway column not found: {pathway_col}")

    if group is not None and cell_type is None:
        current_group_filter = group if isinstance(group, str) else group[0]

        types_to_plot = sorted(data.obs['cell_type'].unique())

        plot_df = data.obs[data.obs['group'] == current_group_filter].copy()
        if plot_df.empty:
            raise ValueError(f"No data in group '{current_group_filter}'")

        types_to_plot = [t for t in types_to_plot if t in plot_df['cell_type'].unique()]
        if not types_to_plot:
             raise ValueError(f"No cell types found in group '{current_group_filter}' to plot.")

        group_scores = []
        for t in types_to_plot:
            scores = plot_df[plot_df['cell_type'] == t][pathway_col].dropna()
            group_scores.append(scores)

        valid_groups = [g for g in group_scores if len(g) >= 3]
        if len(valid_groups) >= 2:
            _, p_global = kruskal(*valid_groups)
        else:
            p_global = 1.0

        all_cell_types = sorted(data.obs['cell_type'].unique())
        type_colors_map = {ct: GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
                           for i, ct in enumerate(all_cell_types)}
        type_palette = {t: type_colors_map[t] for t in types_to_plot if t in type_colors_map}

        fig_width = max(6, len(types_to_plot) * 1.0)
        fig, ax = plt.subplots(figsize=(fig_width, 6))

        sns.violinplot(x='cell_type', y=pathway_col, data=plot_df,
                       order=types_to_plot, palette=type_palette, ax=ax,
                       cut=0, inner=None, linewidth=1.5)

        for i, t in enumerate(types_to_plot):
            median_val = plot_df[plot_df['cell_type'] == t][pathway_col].median()
            if not np.isnan(median_val):
                ax.plot(i, median_val, marker='s', markersize=8, color='white',
                        markeredgecolor='black', markeredgewidth=1.5, zorder=3)

        ax.set_xlabel('')
        ax.set_ylabel('Pathway Score', fontsize=10)
        ax.set_xticklabels(types_to_plot, rotation=45, ha='right', fontsize=8)
        title = f'{simplify_name(pathway_col, method)} in {current_group_filter}'

        legend_elements = [mpatches.Rectangle((0, 0), 1, 1,
                                              facecolor=type_palette[t],
                                              label=t) for t in types_to_plot]
        ax.legend(handles=legend_elements, loc='upper right', frameon=False,
                  fontsize=8, ncol=max(1, len(types_to_plot) // 10))

    elif cell_type is not None and group is None:
        current_type_filter = cell_type if isinstance(cell_type, str) else cell_type[0]

        groups_to_plot = sorted(data.obs['group'].unique())

        plot_df = data.obs[data.obs['cell_type'] == current_type_filter].copy()
        if plot_df.empty:
            raise ValueError(f"No data in cell type '{current_type_filter}'")

        groups_to_plot = [g for g in groups_to_plot if g in plot_df['group'].unique()]
        if not groups_to_plot:
            raise ValueError(f"No groups found for cell type '{current_type_filter}' to plot.")

        group_scores = []
        for g in groups_to_plot:
            scores = plot_df[plot_df['group'] == g][pathway_col].dropna()
            group_scores.append(scores)

        valid_groups = [g for g in group_scores if len(g) >= 3]
        if len(valid_groups) >= 2:
            _, p_global = kruskal(*valid_groups)
        else:
            p_global = 1.0

        all_groups = sorted(data.obs['group'].unique())
        group_colors_map = {g: GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
                            for i, g in enumerate(all_groups)}
        group_palette = {g: group_colors_map[g] for g in groups_to_plot if g in group_colors_map}

        fig, ax = plt.subplots(figsize=(max(4, len(groups_to_plot)*0.8), 6))

        sns.violinplot(x='group', y=pathway_col, data=plot_df, order=groups_to_plot,
                       palette=group_palette, ax=ax, cut=0, inner=None, linewidth=1.5)

        for i, g in enumerate(groups_to_plot):
            median_val = plot_df[plot_df['group'] == g][pathway_col].median()
            if not np.isnan(median_val):
                ax.plot(i, median_val, marker='s', markersize=8, color='white',
                        markeredgecolor='black', markeredgewidth=1.5, zorder=3)

        ax.set_xlabel('')
        ax.set_ylabel('Pathway Score', fontsize=10)
        ax.set_xticklabels(groups_to_plot, rotation=0, fontsize=10)
        title = f'{simplify_name(pathway_col, method)} in {current_type_filter}'

        legend_elements = [mpatches.Rectangle((0, 0), 1, 1,
                                              facecolor=group_palette.get(g, '#CCCCCC'),
                                              label=g) for g in groups_to_plot]
        ax.legend(handles=legend_elements, loc='upper left', frameon=False, fontsize=8)

    else:
        raise ValueError("Please specify either 'group' (to compare types) or 'cell_type' (to compare groups), but not both.")

    y_max = plot_df[pathway_col].max()
    y_min = plot_df[pathway_col].min()
    y_range = y_max - y_min
    if y_range == 0:
        y_range = 1.0

    ax.set_ylim(top=y_max + 0.25 * y_range)

    ax.text(0.5, 0.95, format_pval(p_global), transform=ax.transAxes,
            ha='center', va='top', fontsize=9, color='black')

    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    sns.despine(ax=ax)
    ax.grid(False)

    plt.tight_layout()

    return fig
