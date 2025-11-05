"""
scGeneset/analyze/comm.py - Cell Communication Analysis
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge, FancyArrowPatch
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
import seaborn as sns
import scanpy as sc
from scipy.sparse import issparse
from scipy.interpolate import splprep, splev
from collections import defaultdict
from ..utils.config import CMAP_POS, GROUP_COLORS, GLOBAL_PALETTE
from ..utils.helper import simplify_name


def analyze_communication(data, method='AUCell', group=None, pathway=None, z=True):
   """
   Cell-cell communication analysis

   Parameters:
       data: AnnData object
       method: Scoring method
       group: Group to analyze (default: all cells)
       pathway: If specified, group analysis by high/low scores of this pathway
       z: Whether to use normalized data

   Returns:
       results: Dictionary containing results for each analysis group
           - 'liana_df': LIANA results DataFrame
           - 'fig_heatmap': Heatmap Figure
           - 'fig_chord': Chord diagram Figure
           - 'interaction_matrix': Interaction matrix DataFrame
   """

   try:
       import liana as li
   except ImportError:
       raise ImportError("LIANA required: pip install liana")

   results = {}

   if pathway is not None:
       suffix = '_Z' if z else ''
       pathway_col = f'{method}_{pathway.replace("/", "_")}{suffix}'

       if pathway_col not in data.obs.columns:
           raise ValueError(f"Pathway column not found: {pathway_col}")

       scores = data.obs[pathway_col].dropna()
       if len(scores) < 10:
           raise ValueError("Insufficient pathway score data")

       median_score = scores.median()
       high_cells = data.obs[pathway_col] >= median_score
       low_cells = data.obs[pathway_col] < median_score

       groups_to_analyze = [
           ('High', data[high_cells].copy()),
           ('Low', data[low_cells].copy())
       ]
       print(f"Grouped by {pathway}: High={groups_to_analyze[0][1].n_obs} cells, "
             f"Low={groups_to_analyze[1][1].n_obs} cells")
   else:
       if group is not None:
           groups_to_analyze = [
               (group, data[data.obs['group'] == group].copy())
           ]
       else:
           groups_to_analyze = [('All', data.copy())]

   for group_name, adata_group in groups_to_analyze:
       print(f"\n{'=' * 60}")
       print(f"Analyzing group: {group_name}")
       print(f"{'=' * 60}\n")

       adata_group = _prepare_for_liana(adata_group)

       if adata_group.n_obs == 0:
           print(f"Group {group_name} has no valid cells, skipping")
           continue

       print("Running LIANA analysis...")
       try:
           li.mt.rank_aggregate(
               adata_group,
               groupby='cell_type',
               resource_name='consensus',
               expr_prop=0.05,
               min_cells=5,
               return_all_lrs=False,
               verbose=True,
               use_raw=False,
               n_perms=100
           )
       except Exception as e:
           print(f"LIANA analysis failed: {e}")
           continue

       liana_res = adata_group.uns['liana_res'].copy()

       sig_res = liana_res[liana_res['magnitude_rank'] < 0.05].copy()
       if len(sig_res) == 0:
           sig_res = liana_res[liana_res['magnitude_rank'] < 0.1].copy()
           if len(sig_res) == 0:
               print(f"Group {group_name} has no significant interactions")
               continue

       print(f"Found {len(sig_res)} significant interactions")

       weight_col = None
       for col in ['lr_means', 'lr_weight', 'combined_score', 'score']:
           if col in liana_res.columns:
               weight_col = col
               break
       if weight_col is None:
           weight_col = liana_res.columns[-1]

       fig_heatmap, interaction_matrix = _plot_comm_heatmap(
           sig_res, adata_group, weight_col, group_name
       )

       fig_chord = _plot_chord_diagram(
           sig_res, adata_group, weight_col, group_name
       )

       results[group_name] = {
           'liana_df': sig_res,
           'fig_heatmap': fig_heatmap,
           'fig_chord': fig_chord,
           'interaction_matrix': interaction_matrix
       }

   return results


def _prepare_for_liana(adata):
   """Prepare data for LIANA analysis"""
   adata = adata.copy()
   adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)

   cell_type_counts = adata.obs['cell_type'].value_counts()
   valid_types = cell_type_counts[cell_type_counts >= 5].index.tolist()
   adata = adata[adata.obs['cell_type'].isin(valid_types)].copy()

   if adata.n_obs == 0:
       return adata

   if hasattr(adata, 'raw') and adata.raw is not None:
       raw_adata = adata.raw.to_adata()
       common_genes = adata.var_names.intersection(raw_adata.var_names)
       if len(common_genes) >= adata.n_vars * 0.8:
           raw_subset = raw_adata[:, common_genes].copy()
           adata = adata[:, common_genes].copy()
           adata.X = raw_subset.X.toarray() if issparse(raw_subset.X) else raw_subset.X
       else:
           if issparse(adata.X):
               adata.X = adata.X.toarray()
   else:
       if issparse(adata.X):
           adata.X = adata.X.toarray()

   adata.X = np.nan_to_num(adata.X.astype(float), nan=0.0, posinf=0.0, neginf=0.0)
   adata.X = np.clip(adata.X, 0, None)

   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)
   adata.X = np.nan_to_num(adata.X.astype(float), nan=0.0, posinf=0.0, neginf=0.0)

   min_cells = max(3, int(adata.n_obs * 0.01))
   sc.pp.filter_genes(adata, min_cells=min_cells)

   return adata


def _plot_comm_heatmap(sig_res, adata, weight_col, group_name):
   """Plot communication heatmap"""
   cell_types = sorted(adata.obs['cell_type'].unique())
   interaction_matrix = pd.DataFrame(0.0, index=cell_types, columns=cell_types)

   for _, row in sig_res.iterrows():
       source = row['source']
       target = row['target']
       if source in cell_types and target in cell_types:
           interaction_matrix.loc[source, target] += row[weight_col]

   max_val = interaction_matrix.max().max()
   if max_val > 0:
       interaction_matrix = interaction_matrix / max_val

   fig, ax = plt.subplots(figsize=(10, 8))
   mask = interaction_matrix == 0

   sns.heatmap(interaction_matrix, cmap=CMAP_POS, ax=ax, square=True, mask=mask,
               cbar_kws={'label': 'Normalized Interaction Strength', 'shrink': 0.8},
               linewidths=0, linecolor='white', vmin=0, vmax=1, annot=False)

   ax.set_title(f'Cell-Cell Communication Heatmap\n({group_name})',
                fontsize=12, pad=10)
   ax.set_xlabel('Target Cell Type', fontsize=10)
   ax.set_ylabel('Source Cell Type', fontsize=10)
   plt.xticks(rotation=45, ha='right', fontsize=8)
   plt.yticks(rotation=0, fontsize=8)
   plt.tight_layout()

   return fig, interaction_matrix


def _plot_chord_diagram(sig_res, adata, weight_col, group_name):
   """Plot chord diagram"""
   inflammatory_keywords = {
       'TNF', 'IL1', 'IL6', 'IL8', 'IL17', 'IL23', 'IL12', 'IFNG', 'CXCL',
       'CCL', 'IL18', 'LTA', 'TGFB', 'TGF', 'GMCSF', 'NFKB', 'STAT',
       'RELA', 'JUN', 'FOS', 'MYD88', 'IRF', 'CASP'
   }

   cell_types = sorted(adata.obs['cell_type'].unique())

   send_counts = {}
   receive_counts = {}
   for ct in cell_types:
       send_counts[ct] = len(sig_res[sig_res['source'] == ct])
       receive_counts[ct] = len(sig_res[sig_res['target'] == ct])

   total_counts = {ct: send_counts[ct] + receive_counts[ct] for ct in cell_types}
   total_comm = sum(total_counts.values())

   if total_comm == 0:
       print("No communication activity, cannot plot chord diagram")
       return None

   total_chord_angle = 2 * np.pi * 0.75
   angles = {ct: (total_counts[ct] / total_comm) * total_chord_angle
             for ct in cell_types}

   cell_angles = {}
   current_angle = (2 * np.pi - total_chord_angle) / 2

   for ct in cell_types:
       start = current_angle
       end = current_angle + angles[ct]
       cell_angles[ct] = {
           'start': start,
           'end': end,
           'mid': (start + end) / 2,
           'send': send_counts[ct],
           'receive': receive_counts[ct]
       }
       current_angle = end + (2 * np.pi * 0.25 / len(cell_types))

   fig, ax = plt.subplots(figsize=(16, 16))
   ax.set_xlim(-2.0, 2.0)
   ax.set_ylim(-2.0, 2.0)
   ax.set_aspect('equal')
   ax.axis('off')

   radius = 1.0
   arc_width = 0.12

   def get_arc_color(ct):
       send = cell_angles[ct]['send']
       receive = cell_angles[ct]['receive']
       total = send + receive
       if total == 0:
           return '#CCCCCC'
       send_ratio = send / total
       if send_ratio > 0.6:
           return '#C2ABC8'
       elif send_ratio < 0.4:
           return '#A3C9D5'
       else:
           return '#F2BB6B'

   for ct in cell_types:
       angles_data = cell_angles[ct]
       start_deg = np.degrees(angles_data['start'])
       end_deg = np.degrees(angles_data['end'])
       color = get_arc_color(ct)

       wedge = Wedge((0, 0), radius + arc_width, start_deg, end_deg,
                     width=arc_width, facecolor=color, edgecolor='black',
                     linewidth=2.5, zorder=3)
       ax.add_patch(wedge)

       mid_angle = angles_data['mid']
       label_radius = radius + arc_width / 2 + 0.25
       label_x = label_radius * np.cos(mid_angle)
       label_y = label_radius * np.sin(mid_angle)

       rotation = np.degrees(mid_angle)
       if rotation > 90 and rotation < 270:
           rotation = rotation + 180
           ha = 'right'
       else:
           ha = 'left'

       label_text = f"{ct}\n(S:{angles_data['send']}, R:{angles_data['receive']})"
       ax.text(label_x, label_y, label_text, ha=ha, va='center',
               fontsize=9, fontweight='bold', rotation=rotation,
               rotation_mode='anchor')

   connection_data = []

   ligand_col = None
   receptor_col = None
   for col in sig_res.columns:
       if 'ligand' in col.lower():
           ligand_col = col
       if 'receptor' in col.lower():
           receptor_col = col

   for _, row in sig_res.iterrows():
       source = row['source']
       target = row['target']
       if source not in cell_types or target not in cell_types:
           continue

       is_inflammatory = False
       if ligand_col and pd.notna(row[ligand_col]):
           if any(kw in str(row[ligand_col]).upper() for kw in inflammatory_keywords):
               is_inflammatory = True
       if receptor_col and pd.notna(row[receptor_col]) and not is_inflammatory:
           if any(kw in str(row[receptor_col]).upper() for kw in inflammatory_keywords):
               is_inflammatory = True

       connection_data.append({
           'source': source,
           'target': target,
           'weight': row[weight_col],
           'inflammatory': is_inflammatory
       })

   connection_data = sorted(connection_data, key=lambda x: x['weight'])

   outgoing_points = defaultdict(list)
   incoming_points = defaultdict(list)

   for conn in connection_data:
       outgoing_points[conn['source']].append(conn)
       incoming_points[conn['target']].append(conn)

   inner_radius = radius

   for ct in cell_types:
       outs = outgoing_points[ct]
       if outs:
           sorted_outs = sorted(outs, key=lambda c: cell_angles[c['target']]['mid'])
           angle_points = np.linspace(cell_angles[ct]['start'],
                                      cell_angles[ct]['end'],
                                      len(outs) + 2)[1:-1]
           for idx, conn in enumerate(sorted_outs):
               conn['source_angle'] = angle_points[idx]

       ins = incoming_points[ct]
       if ins:
           sorted_ins = sorted(ins, key=lambda c: cell_angles[c['source']]['mid'])
           angle_points = np.linspace(cell_angles[ct]['start'],
                                      cell_angles[ct]['end'],
                                      len(ins) + 2)[1:-1]
           for idx, conn in enumerate(sorted_ins):
               conn['target_angle'] = angle_points[idx]

   for conn in connection_data:
       source_angle = conn.get('source_angle', cell_angles[conn['source']]['mid'])
       target_angle = conn.get('target_angle', cell_angles[conn['target']]['mid'])

       source_x = inner_radius * np.cos(source_angle)
       source_y = inner_radius * np.sin(source_angle)
       target_x = inner_radius * np.cos(target_angle)
       target_y = inner_radius * np.sin(target_angle)

       angle_diff = abs(source_angle - target_angle)
       if angle_diff > np.pi:
           angle_diff = 2 * np.pi - angle_diff

       mid_angle = (source_angle + target_angle) / 2
       if abs(target_angle - source_angle) > np.pi:
           mid_angle += np.pi

       depth_factor = 0.3 + 0.4 * (angle_diff / np.pi)
       control_radius = inner_radius * (1 - depth_factor)
       control_x = control_radius * np.cos(mid_angle)
       control_y = control_radius * np.sin(mid_angle)

       points = np.array([[source_x, source_y], [control_x, control_y],
                          [target_x, target_y]])
       tck, u = splprep(points.T, s=0, k=2)
       u_new = np.linspace(0, 1, 100)
       x_new, y_new = splev(u_new, tck)

       segments = [np.column_stack([x_new[i:i + 2], y_new[i:i + 2]])
                   for i in range(len(x_new) - 1)]

       if conn['inflammatory']:
           cmap = mcolors.LinearSegmentedColormap.from_list(
               'inflammatory_cmap', ["#AAD09D", "#FFF7AC", "#ECB477"]
           )
           alpha = 0.8
           zorder = 2
       else:
           cmap = mcolors.LinearSegmentedColormap.from_list(
               'general_cmap', ["#D5E8D0", "#FFF9D9"]
           )
           alpha = 0.3
           zorder = 1

       lc = LineCollection(segments, cmap=cmap,
                           norm=mcolors.Normalize(0, 1),
                           linewidth=conn['weight'] * 4, alpha=alpha,
                           zorder=zorder)
       lc.set_array(np.linspace(0, 1, len(segments)))
       ax.add_collection(lc)

   from matplotlib.lines import Line2D

   legend_elements = [
       mpatches.Rectangle((0, 0), 1, 1, facecolor='#C2ABC8',
                          edgecolor='black', linewidth=1.5,
                          label='Sender (>60% sent)'),
       mpatches.Rectangle((0, 0), 1, 1, facecolor='#A3C9D5',
                          edgecolor='black', linewidth=1.5,
                          label='Receiver (<40% sent)'),
       mpatches.Rectangle((0, 0), 1, 1, facecolor='#F2BB6B',
                          edgecolor='black', linewidth=1.5,
                          label='Balanced (40-60%)'),
       mpatches.Rectangle((0, 0), 1, 1, facecolor='none',
                          edgecolor='none', label=''),
       Line2D([0], [0], color='#D5E8D0', lw=3, alpha=0.5,
              label='General'),
       Line2D([0], [0], color='#ECB477', lw=4, alpha=0.8,
              label='Inflammatory')
   ]

   ax.legend(handles=legend_elements, loc='center right',
             bbox_to_anchor=(1.2, 0.5), fontsize=10, frameon=True,
             fancybox=True, shadow=True, framealpha=0.95)

   cbar_ax = ax.inset_axes([0.85, 0.05, 0.03, 0.15])
   cbar_cmap = mcolors.LinearSegmentedColormap.from_list(
       'chord_gradient', ["#AAD09D", "#FFF7AC", "#ECB477"]
   )
   norm = mcolors.Normalize(vmin=0, vmax=1)
   cb = plt.colorbar(
       plt.cm.ScalarMappable(norm=norm, cmap=cbar_cmap),
       cax=cbar_ax, orientation='vertical'
   )
   cb.set_label('Interaction\nStrength', rotation=0, ha='left', va='center',
                fontsize=8, labelpad=10)
   cb.ax.tick_params(labelsize=7)
   cb.set_ticks([0, 0.5, 1])
   cb.set_ticklabels(['Weak', 'Medium', 'Strong'])

   title = 'Cell-Cell Communication Chord Diagram\n(Arc length proportional to Total communication activity)'
   ax.text(0, -1.75, title, ha='center', va='top',
           fontsize=13, fontweight='bold')

   explanation = 'S: Ligands sent | R: Ligands received | Arcs bent inward toward center'
   ax.text(0, -1.88, explanation, ha='center', va='top',
           fontsize=9, style='italic', color='gray')

   plt.tight_layout()

   return fig
