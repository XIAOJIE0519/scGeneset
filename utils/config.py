"""
scInflammation/utils/config.py - Global configuration and constants
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ========== Global plotting settings ==========
def setup_plotting():
    """Set global plotting parameters"""
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['legend.fontsize'] = 8
    plt.rcParams['legend.title_fontsize'] = 10
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['lines.linewidth'] = 1.5
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['savefig.dpi'] = 600
    plt.rcParams['savefig.transparent'] = False


# ========== Color schemes ==========
def create_global_palette():
    """Create a palette of 100 colors"""
    palette_tab20 = sns.color_palette("tab20", 20)
    palette_tab20b = sns.color_palette("tab20b", 20)
    palette_tab20c = sns.color_palette("tab20c", 20)
    palette_set3 = sns.color_palette("Set3", 12)
    palette_paired = sns.color_palette("Paired", 12)

    combined = (list(palette_tab20) + list(palette_tab20b) +
                list(palette_tab20c) + list(palette_set3) + list(palette_paired))

    if len(combined) < 100:
        remaining = 100 - len(combined)
        palette_hls = sns.color_palette("hls", n_colors=remaining)
        combined.extend(palette_hls)

    return combined[:100]


GLOBAL_PALETTE = create_global_palette()
GROUP_COLORS = {}
CMAP_NEG_POS = 'RdBu_r'
CMAP_POS = mcolors.LinearSegmentedColormap.from_list(
    'custom_pos', ["#AAD09D", "#ECF4DD", "#FFF7AC", "#ECB477"]
)
CMAP_GREEN_RED = mcolors.LinearSegmentedColormap.from_list(
    'green_red', ["green", "yellow", "red"]
)

LEGEND_KWARGS = {
    'frameon': False, 'fontsize': 8, 'title_fontsize': 10,
    'markerscale': 1.0, 'handlelength': 1.5, 'handletextpad': 0.5
}

# ========== Parallel computing ==========
N_JOBS = -1

# ========== Pathway configuration ==========
PATHWAY_FILE = 'pathway.csv'  # Will be placed in the library's data directory

setup_plotting()
