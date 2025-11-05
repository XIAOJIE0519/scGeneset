"""
scGeneset/utils/__init__.py
"""

from .config import (
    setup_plotting,
    GLOBAL_PALETTE,
    GROUP_COLORS,
    CMAP_POS,
    CMAP_NEG_POS,
    CMAP_GREEN_RED,
    LEGEND_KWARGS,
    N_JOBS,
    PATHWAY_FILE,
)

from .helper import (
    simplify_name,
    format_pval,
    find_circle,
    calc_gene_score_subset,
    get_residuals,
)

__all__ = [
    'setup_plotting',
    'GLOBAL_PALETTE',
    'GROUP_COLORS',
    'CMAP_POS',
    'CMAP_NEG_POS',
    'CMAP_GREEN_RED',
    'LEGEND_KWARGS',
    'N_JOBS',
    'PATHWAY_FILE',
    'simplify_name',
    'format_pval',
    'find_circle',
    'calc_gene_score_subset',
    'get_residuals',
]
