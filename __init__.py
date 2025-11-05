"""
scGeneset: A toolkit for single-cell gene set analysis

Main features:
- Data loading and preprocessing
- Multiple pathway scoring methods (AUCell, AddModule, ssGSEA, singscore)
- Visualization (heatmap, violin plot, UMAP, etc.)
- Advanced analysis (consistency, WGCNA, cell communication, pseudotime, etc.)
"""

__version__ = '0.0.1'
__author__ = 'Shanjie Luan'
__email__ = 'Luan20050519@163.com'

# ===========================
# 1. Import core functions
# ===========================
from .core.loader import load_data, load_pathway_genes
from .core.scorer import score_pathways

# ===========================
# 2. Import plotting functions
# ===========================
from .plot.heat import plot_heatmap
from .plot.violin import plot_violin
from .plot.umap import plot_umap
from .plot.corr import plot_corr_scatter, plot_corr_matrix
from .plot.pca import plot_pca

# ===========================
# 3. Import analysis functions
# ===========================
from .analyze.consist import calc_consistency
from .analyze.wgcna import calc_wgcna
from .analyze.comm import analyze_communication
from .analyze.pseudo import analyze_pseudotime
from .analyze.reclust import recluster_cells
from .analyze.diff import analyze_differential

# ===========================
# 4. Import utility functions
# ===========================
from .utils.config import (
    setup_plotting,
    GLOBAL_PALETTE,
    GROUP_COLORS,
    N_JOBS,
    PATHWAY_FILE,
)

from .utils.helper import (
    simplify_name,
    format_pval,
    calc_gene_score_subset,
)

# ===========================
# Define public API
# ===========================
__all__ = [
    # Core functions
    'load_data',
    'load_pathway_genes',
    'score_pathways',

    # Plotting functions
    'plot_heatmap',
    'plot_violin',
    'plot_umap',
    'plot_corr_scatter',
    'plot_corr_matrix',
    'plot_pca',

    # Analysis functions
    'calc_consistency',
    'calc_wgcna',
    'analyze_communication',
    'analyze_pseudotime',
    'recluster_cells',
    'analyze_differential',

    # Utility functions
    'setup_plotting',
    'simplify_name',
    'format_pval',
    'calc_gene_score_subset',

    # Configuration constants
    'GLOBAL_PALETTE',
    'GROUP_COLORS',
    'PATHWAY_FILE',
]

setup_plotting()

print("✓ scInflammation imported successfully!")
print(f"Library location: {__file__}")

# Display current working directory for debugging
import os
print(f"Current working directory: {os.getcwd()}")
