"""
scInflammation/utils/helper.py
"""
import numpy as np
import pandas as pd
import math
from scipy.stats import spearmanr
from joblib import Parallel, delayed
from .config import N_JOBS


def simplify_name(pathway, method=''):
    """Simplify pathway name"""
    pathway = pathway.replace(f'{method}_', '')
    return (pathway.replace('_R_HSA_168256', '')
            .replace('_GO0002685', '')
            .replace('_GO0002697', '')
            .replace('_', ' ').title())


def format_pval(p):
    """Format p-value"""
    if p < 0.001:
        return "p < 0.001"
    elif p < 0.01:
        return f"p = {p:.3f}"
    else:
        return f"p = {p:.2f}"


def find_circle(p1, p2, p3):
    """Find circle through three points"""
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    a = 2 * (x2 - x1)
    b = 2 * (y2 - y1)
    c = x2 ** 2 + y2 ** 2 - x1 ** 2 - y1 ** 2
    d = 2 * (x3 - x1)
    e = 2 * (y3 - y1)
    f = x3 ** 2 + y3 ** 2 - x1 ** 2 - y1 ** 2
    determinant = a * e - b * d
    if abs(determinant) < 1e-6:
        return None
    cx = (c * e - b * f) / determinant
    cy = (a * f - d * c) / determinant
    radius = math.sqrt((cx - x1) ** 2 + (cy - y1) ** 2)
    return (cx, cy), radius


def calc_gene_score_subset(adata_subset, gene_set, method='AUCell'):
    """Calculate gene set score for subset"""
    from .config import N_JOBS

    mat = adata_subset.raw.to_adata()
    gene_exp_array = mat.X.toarray() if hasattr(mat.X, 'toarray') else mat.X
    gene_names = mat.var_names.tolist()
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    indices = [gene_to_idx[g] for g in gene_set if g in gene_to_idx]
    if len(indices) < 2:
        return np.zeros(adata_subset.n_obs)

    indices_np = np.array(indices)
    n_cells = adata_subset.n_obs
    n_genes = mat.n_vars

    if method == 'AUCell':
        from ..core.scorer import aucell_score_single_cell
        scores = Parallel(n_jobs=N_JOBS, backend='threading')(
            delayed(aucell_score_single_cell)(i, gene_exp_array, indices_np)
            for i in range(n_cells)
        )
    elif method == 'AddModule':
        from ..core.scorer import seurat_score_single
        _, scores = seurat_score_single("gene_set", list(gene_set), mat)
    elif method == 'ssGSEA':
        from ..core.scorer import ssgsea_score_single_cell
        scores = Parallel(n_jobs=N_JOBS, backend='threading')(
            delayed(ssgsea_score_single_cell)(i, gene_exp_array, indices_np, n_genes)
            for i in range(n_cells)
        )
    elif method == 'singscore':
        from ..core.scorer import calc_singscore_optimized
        scores = Parallel(n_jobs=N_JOBS, backend='threading')(
            delayed(calc_singscore_optimized)(i, gene_exp_array, indices_np, n_genes)
            for i in range(n_cells)
        )
    else:
        scores = np.zeros(n_cells)

    return np.array(scores)


def get_residuals(score_vector, covariate_vector):
    """Calculate OLS regression residuals"""
    import statsmodels.api as sm
    df = pd.DataFrame({'Y': score_vector, 'X': covariate_vector}).dropna()
    if df.shape[0] < 10 or df['X'].std() == 0:
        return pd.Series(np.nan, index=score_vector.index)

    X = sm.add_constant(df['X'])
    Y = df['Y']
    model = sm.OLS(Y, X).fit()
    residuals = pd.Series(model.resid, index=df.index)
    return residuals.reindex(score_vector.index)
