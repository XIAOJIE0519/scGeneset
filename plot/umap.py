"""
scGeneset/plot/umap.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import seaborn as sns
import scanpy as sc
from ..utils.config import CMAP_POS, GROUP_COLORS, GLOBAL_PALETTE
from ..utils.helper import simplify_name


def plot_umap(data, method='AUCell', z=True, top_n=5):
    """
    Plot UMAP with circular pathway annotations.

    Args:
        data: AnnData object
        method: Scoring method
        z: Whether to use Z-normalized data
        top_n: Number of top pathways to display in inner ring

    Returns:
        fig: matplotlib Figure object
    """

    if 'X_umap' not in data.obsm:
        sc.tl.umap(data)

    umap_coords = data.obsm['X_umap']

    cell_types = sorted([ct for ct in data.obs['cell_type'].unique()
                         if data.obs[data.obs['cell_type'] == ct].shape[0] > 0])
    num_cell_types = len(cell_types)

    suffix = '_Z' if z else ''
    pathway_cols = [col for col in data.obs.columns
                    if col.startswith(f'{method}_') and col.endswith(suffix)
                    and 'Custom' not in col][:top_n]

    if len(pathway_cols) == 0:
        raise ValueError(f"No valid pathways found (method={method}, z={z})")

    cell_type_categories = sorted(data.obs['cell_type'].unique())
    cell_type_colors = {ct: GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
                        for i, ct in enumerate(cell_type_categories)}

    group_categories = sorted(data.obs['group'].unique())
    group_palette = {
        g: GLOBAL_PALETTE[i % len(GLOBAL_PALETTE)]
        for i, g in enumerate(group_categories)
    }

    cell_type_data = {}
    total_cell_counts = {}

    for ct in cell_types:
        subset = data[data.obs['cell_type'] == ct]
        group_counts = subset.obs['group'].value_counts()
        total_cells = group_counts.sum()
        total_cell_counts[ct] = total_cells

        group_props = {}
        for g in group_categories:
            group_props[g] = group_counts.get(g, 0) / total_cells if total_cells > 0 else 0

        pathway_scores_avg = {}
        for p_col in pathway_cols:
            if p_col in subset.obs.columns:
                pathway_scores_avg[p_col] = subset.obs[p_col].mean()
            else:
                pathway_scores_avg[p_col] = 0

        cell_type_data[ct] = {
            'color': cell_type_colors[ct],
            'group_proportions': group_props,
            'pathway_scores_avg': pathway_scores_avg,
        }

    all_scores = []
    for ct_data in cell_type_data.values():
        all_scores.extend(ct_data['pathway_scores_avg'].values())

    min_score = np.min(all_scores) if all_scores else 0
    max_score = np.max(all_scores) if all_scores else 1
    norm_pathway = mcolors.Normalize(vmin=min_score, vmax=max_score)

    fig = plt.figure(figsize=(18, 16))
    fig.patch.set_facecolor('white')

    fig_width, fig_height = fig.get_size_inches()

    ax_size_inches = 10
    ax_width_fig = ax_size_inches / fig_width
    ax_height_fig = ax_size_inches / fig_height
    ax_left_margin = 0.05
    ax_bottom_margin = (1 - ax_height_fig) / 2

    ax_umap_bbox = [ax_left_margin, ax_bottom_margin, ax_width_fig, ax_height_fig]
    ax_umap = fig.add_axes(ax_umap_bbox)
    ax_umap.set_facecolor('white')

    umap_x_min, umap_x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
    umap_y_min, umap_y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()
    umap_x_range = umap_x_max - umap_x_min
    umap_y_range = umap_y_max - umap_y_min
    max_range = max(umap_x_range, umap_y_range)
    umap_x_center = (umap_x_min + umap_x_max) / 2
    umap_y_center = (umap_y_min + umap_y_max) / 2

    padding_factor = 1.1
    data_radius = 0.5 / padding_factor
    ring_width = 0.02
    gap_width = 0.01

    radius_inner = data_radius + gap_width + ring_width
    radius_middle = radius_inner + gap_width + ring_width
    radius_outer = radius_middle + gap_width + ring_width

    umap_radius = (radius_inner - ring_width) * 1.414

    umap_coords_transformed = np.zeros_like(umap_coords)
    umap_coords_transformed[:, 0] = (umap_coords[:, 0] - umap_x_center) / max_range * umap_radius + 0.5
    umap_coords_transformed[:, 1] = (umap_coords[:, 1] - umap_y_center) / max_range * umap_radius + 0.5

    ax_umap.set_xlim(0, 1)
    ax_umap.set_ylim(0, 1)
    ax_umap.set_aspect('equal')

    for ct in cell_types:
        mask = data.obs['cell_type'] == ct
        ax_umap.scatter(umap_coords_transformed[mask, 0], umap_coords_transformed[mask, 1],
                        color=cell_type_colors[ct], s=20, alpha=0.3,
                        label=ct, edgecolors='none',
                        transform=ax_umap.transAxes)

    sns.kdeplot(x=umap_coords_transformed[:, 0], y=umap_coords_transformed[:, 1],
                ax=ax_umap, color='gray', levels=5, linewidths=1.5,
                alpha=0.7, linestyles='dashed',
                transform=ax_umap.transAxes)

    ax_umap.set_xlabel('UMAP_1', fontsize=10)
    ax_umap.set_ylabel('UMAP_2', fontsize=10)
    ax_umap.set_xticks([])
    ax_umap.set_yticks([])
    ax_umap.spines['top'].set_visible(False)
    ax_umap.spines['right'].set_visible(False)
    ax_umap.spines['bottom'].set_visible(False)
    ax_umap.spines['left'].set_visible(False)
    z_label = " (Z-normalized)" if z else ""
    ax_umap.set_title(f'UMAP with Pathway Scores\n[{method}{z_label}]',
                      fontsize=12, pad=10)

    center_x, center_y = (0.5, 0.5)
    segment_gap = 0.072
    total_angle = 2 * np.pi
    angle_per_segment = total_angle / num_cell_types
    current_angle = 0

    for i, ct in enumerate(cell_types):
        data_ct = cell_type_data[ct]
        start_angle = current_angle + segment_gap / 2
        end_angle = current_angle + angle_per_segment - segment_gap / 2
        start_deg = np.degrees(start_angle)
        end_deg = np.degrees(end_angle)

        wedge_outer = mpatches.Wedge((center_x, center_y), radius_outer,
                                     start_deg, end_deg, width=ring_width,
                                     facecolor=data_ct['color'], edgecolor='black',
                                     linewidth=1.0, zorder=2,
                                     transform=ax_umap.transAxes, clip_on=False)
        ax_umap.add_artist(wedge_outer)

        current_group_angle = start_angle
        for group_name in group_categories:
            prop = data_ct['group_proportions'].get(group_name, 0)
            group_end = current_group_angle + prop * (end_angle - start_angle)
            if prop > 0:
                wedge_mid = mpatches.Wedge((center_x, center_y), radius_middle,
                                           np.degrees(current_group_angle),
                                           np.degrees(group_end), width=ring_width,
                                           facecolor=group_palette[group_name],
                                           edgecolor='black', linewidth=1.0, zorder=2,
                                           transform=ax_umap.transAxes, clip_on=False)
                ax_umap.add_artist(wedge_mid)
            current_group_angle = group_end

        if pathway_cols:
            pathway_angle = (end_angle - start_angle) / len(pathway_cols)
            current_pathway_angle = start_angle
            for p_col in pathway_cols:
                score = data_ct['pathway_scores_avg'][p_col]
                color = CMAP_POS(norm_pathway(score))
                pathway_end = current_pathway_angle + pathway_angle
                wedge_inner = mpatches.Wedge((center_x, center_y), radius_inner,
                                             np.degrees(current_pathway_angle),
                                             np.degrees(pathway_end), width=ring_width,
                                             facecolor=color, edgecolor='black',
                                             linewidth=1.0, zorder=2,
                                             transform=ax_umap.transAxes, clip_on=False)
                ax_umap.add_artist(wedge_inner)
                current_pathway_angle = pathway_end

        current_angle += angle_per_segment

    legend_left_margin = ax_left_margin + ax_width_fig + 0.08
    legend_width = 1.0 - legend_left_margin - 0.02
    ax_legend_bbox = [legend_left_margin, ax_bottom_margin, legend_width, ax_height_fig]
    ax_legend = fig.add_axes(ax_legend_bbox)
    ax_legend.axis('off')

    y_pos = 0.95
    y_step = 0.025

    ax_legend.text(0, y_pos, 'Cell Types (Count)', fontsize=10, fontweight='bold',
                   va='top', ha='left', transform=ax_legend.transAxes)
    y_pos -= 0.03

    for idx, ct in enumerate(cell_types):
        count = total_cell_counts[ct]
        ax_legend.plot(0.01, y_pos - idx * y_step, 'o',
                       color=cell_type_colors[ct], markersize=6,
                       transform=ax_legend.transAxes)
        ax_legend.text(0.03, y_pos - idx * y_step, f'{ct} ({count})',
                       fontsize=8, va='center', ha='left', color='black',
                       transform=ax_legend.transAxes)

    y_pos = y_pos - (len(group_categories) + 1) * y_step

    ax_legend.text(0, y_pos, 'Groups', fontsize=10, fontweight='bold',
                   va='top', ha='left', transform=ax_legend.transAxes)
    y_pos -= 0.03

    for idx, g in enumerate(group_categories):
        ax_legend.add_patch(mpatches.Rectangle((0, y_pos - idx * y_step), 0.02, 0.02,
                                               facecolor=group_palette[g],
                                               edgecolor='black', linewidth=0.5,
                                               transform=ax_legend.transAxes))
        ax_legend.text(0.03, y_pos - idx * y_step + 0.01, g, fontsize=8,
                       va='center', ha='left', transform=ax_legend.transAxes)

    ax_legend.text(0, y_pos - (len(group_categories) + 1) * y_step,
                   'Arc length represents proportion', fontsize=7,
                   va='top', ha='left', style='italic', color='gray',
                   transform=ax_legend.transAxes)

    y_pos = y_pos - (len(group_categories) + 2) * y_step

    ax_legend.text(0, y_pos, f'Pathway Scores (Top {len(pathway_cols)})',
                   fontsize=10, fontweight='bold', va='top', ha='left',
                   transform=ax_legend.transAxes)

    cbar_ax = ax_legend.inset_axes([0.05, y_pos - 0.25, 0.03, 0.2])
    sm = plt.cm.ScalarMappable(cmap=CMAP_POS, norm=norm_pathway)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Normalized Score', rotation=270, labelpad=10, fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    pathway_label_y_start = y_pos - 0.05
    for idx, p_col in enumerate(pathway_cols):
        p_name = simplify_name(p_col.replace('_Z', ''), method)
        ax_legend.text(0.15, pathway_label_y_start - idx * y_step,
                       f'- {p_name}', fontsize=7, va='center', ha='left',
                       transform=ax_legend.transAxes)

    plt.tight_layout()

    return fig
