# scGeneset：单细胞基因集分析流程

## 简介

`scGeneset` 是一个用于单细胞 RNA 测序数据基因集分析的 Python 库。它提供了一整套从数据加载、预处理、基因集通路打分到高级分析和可视化的功能。旨在帮助研究人员轻松探索单细胞数据中的通路活性和细胞状态。

**主要功能:**

*   **数据管理:** 支持从 `h5ad` 文件或独立的 CSV/TSV 文件加载和预处理单细胞数据 (`AnnData` 对象)。
*   **基因集打分:** 实现了多种主流的基因集富集打分方法，包括 AUCell, AddModuleScore (Seurat 风格), ssGSEA 和 singscore。
*   **可视化:** 提供丰富的图表类型，如通路活性热图、小提琴图、UMAP、PCA 和相关性矩阵，直观展示分析结果。
*   **高级分析:** 集成了多种高级分析模块，如通路一致性分析、WGCNA (通路共表达网络)、细胞通讯、拟时序分析、重聚类和差异基因分析。

## 输入数据格式（非常重要！）

`scGeneset` 库的核心数据结构是 `anndata.AnnData` 对象。`load_data` 函数是数据输入的起点，它支持两种主要的数据输入方式：

1.  **首选方式：单个 `h5ad` 文件**
    如果你的数据已经是一个 `AnnData` 对象并保存为 `.h5ad` 文件，这是最推荐的方式。
    该 `h5ad` 文件必须包含：
    *   `.X`: 原始或处理后的表达矩阵。
    *   `.obs`: 包含细胞元数据的 DataFrame。至少需要包含：
        *   **细胞类型（Cell Type）信息**: 默认列名为 `'final_annotation'` (可通过 `annotation_col` 参数指定)。
        *   **分组（Group）信息**: 默认列名为 `'group'` (可通过 `group_col` 参数指定)，例如处理组、对照组或不同的样本。
    *   `.var`: 包含基因元数据的 DataFrame。
    *   `.raw`: **强烈推荐**包含原始计数数据作为 `.raw` 属性，`scGeneset` 中的许多打分和分析功能都需要访问原始计数。如果不存在，`load_data` 会尝试创建一个副本。
    *   `.obsm['X_umap']` 或 `.obsm['X_pca']`: 用于某些可视化（如 UMAP、PCA）的降维嵌入。如果不存在，相关函数会尝试计算。

2.  **备选方式：三个独立的 CSV/TSV/H5AD 文件**
    如果你没有现成的 `h5ad` 文件，你可以提供以下三个文件：
    *   **原始表达矩阵文件 (`raw_file`)**:
        *   **格式**: CSV (`.csv`), TSV (`.txt`, `.tsv`) 或 `h5ad`。
        *   **内容**: 行为细胞 (Cell ID)，列为基因 (Gene ID) 的原始计数矩阵。
        *   **索引/列名**: 细胞 ID 和基因 ID 必须与其他文件中的名称完全匹配。
    *   **细胞注释文件 (`annotation_file`)**:
        *   **格式**: CSV (`.csv`), TSV (`.txt`, `.tsv`)。
        *   **内容**: 行为细胞 ID，第一列为每个细胞的细胞类型注释。
        *   **索引**: 细胞 ID 必须与 `raw_file` 的行名匹配。
    *   **分组信息文件 (`group_file`)**:
        *   **格式**: CSV (`.csv`), TSV (`.txt`, `.tsv`)。
        *   **内容**: 行为细胞 ID，第一列为每个细胞所属的实验组/条件。
        *   **索引**: 细胞 ID 必须与 `raw_file` 的行名匹配。

    `load_data` 函数将使用这些文件来构建 `AnnData` 对象，并自动将其 `X` 矩阵设置为处理后的数据，`.raw` 属性设置为原始计数。

### `AnnData` 对象中的关键列约定

为了方便库内函数的使用，`load_data` 会确保 `adata.obs` 中存在以下两列：
*   `adata.obs['cell_type']`: 存储细胞类型信息，通常是 `annotation_col` 指定的列的副本，并转换为 `category` 类型。
*   `adata.obs['group']`: 存储分组信息，通常是 `group_col` 指定的列的副本，并转换为 `category` 类型。
请确保你的输入数据能够提供这些信息。

---

## 核心工作流程概览

一个典型的 `scGeneset` 分析会遵循以下步骤：

1.  **加载数据**: 使用 `load_data` 函数将原始数据或已处理的 `AnnData` 对象加载到内存中。
2.  **加载通路基因集**: 使用 `load_pathway_genes` 函数加载预定义的或自定义的基因集列表。
3.  **计算通路得分**: 使用 `score_pathways` 函数对 `AnnData` 对象中的每个细胞计算每个基因集的通路活性得分。得分会添加到 `adata.obs` 中。
4.  **进行分析与可视化**: 使用 `plot.*` 和 `analyze.*` 模块中的函数生成各种图表和结果。

## 函数参考

### `scGeneset.core.loader`

#### `load_data(annotation_file=None, group_file=None, raw_file=None, adata_file=None, annotation_col='final_annotation', group_col='group', min_genes=None, min_cells=None, target_sum=None, n_top_genes=None, log_normalize=False, verbose=True)`

*   **功能**: 加载单细胞数据，进行基本预处理，并正确设置 `AnnData` 对象的 `.raw` 属性。
*   **参数**:
    *   `annotation_file` (str, 可选): 细胞注释文件路径 (CSV/TSV)。如果提供了 `adata_file` 则忽略。
    *   `group_file` (str, 可选): 细胞分组文件路径 (CSV/TSV)。如果提供了 `adata_file` 则忽略。
    *   `raw_file` (str, 可选): 原始表达矩阵文件路径 (CSV/TSV/H5AD)。如果提供了 `adata_file` 则忽略。
    *   `adata_file` (str, 可选): 现有的 `h5ad` 文件路径。**如果提供，则优先使用，其他文件路径参数将被忽略。**
    *   `annotation_col` (str, 默认值: `'final_annotation'`): `adata.obs` 中用于细胞类型注释的列名。
    *   `group_col` (str, 默认值: `'group'`): `adata.obs` 中用于分组信息的列名。
    *   `min_genes` (int, 可选): 细胞至少表达的基因数。低于此值的细胞将被过滤。
    *   `min_cells` (int, 可选): 基因至少在多少个细胞中表达。低于此值的基因将被过滤。
    *   `target_sum` (float, 可选): 归一化后每个细胞的总计数。例如，`1e4`。
    *   `n_top_genes` (int, 可选): 要选择的高度可变基因的数量。如果为 `None` 则不选择。
    *   `log_normalize` (bool, 默认值: `False`): 在归一化后是否应用 `log1p` 转换。
    *   `verbose` (bool, 默认值: `True`): 是否打印详细的加载和预处理信息。
*   **返回**: `anndata.AnnData` 对象。
*   **示例（h5ad数据可在data文件夹底下找到）**:

    ```python
    import scGeneset as sci
    import os

    # 示例1: 从 h5ad 文件加载 (推荐)
    # 假设 'annotated_adata.h5ad' 包含 .X, .obs (有'final_annotation'和'group'列), .var, .raw
    adata_h5ad = sci.load_data(adata_file='./单细胞输出/初步分析/annotated_adata.h5ad', verbose=True)
    print(f"Loaded from h5ad: {adata_h5ad.n_obs} cells, {adata_h5ad.n_vars} genes.")

    # 示例2: 从独立文件加载 (需要准备好三个文件)
    # 创建一些假数据文件
    # if not os.path.exists('temp_raw.csv'):
    #     import pandas as pd
    #     import numpy as np
    #     cells = [f'cell_{i}' for i in range(100)]
    #     genes = [f'gene_{i}' for i in range(500)]
    #     pd.DataFrame(np.random.randint(0, 100, size=(100, 500)), index=cells, columns=genes).to_csv('temp_raw.csv')
    #     pd.DataFrame({'cell_type': np.random.choice(['T_cell', 'B_cell'], 100)}, index=cells).to_csv('temp_annotation.csv')
    #     pd.DataFrame({'group': np.random.choice(['Ctrl', 'Treat'], 100)}, index=cells).to_csv('temp_group.csv')

    # adata_files = sci.load_data(
    #     annotation_file='temp_annotation.csv',
    #     group_file='temp_group.csv',
    #     raw_file='temp_raw.csv',
    #     annotation_col='cell_type', # 假设文件中列名是'cell_type'
    #     group_col='group',         # 假设文件中列名是'group'
    #     min_genes=200,
    #     target_sum=1e4,
    #     log_normalize=True,
    #     verbose=True
    # )
    # print(f"Loaded from files: {adata_files.n_obs} cells, {adata_files.n_vars} genes.")
    ```
    

### `scGeneset.core.scorer`

#### `score_pathways(adata, pathway_dict=None, custom_genes=None, methods=['AUCell', 'AddModule', 'ssGSEA', 'singscore'], normalize=True, min_genes=2)`

*   **功能**: 计算 `AnnData` 对象中每个细胞的基因集通路活性得分，并将结果添加到 `adata.obs`。
*   **参数**:
    *   `adata` (`AnnData`, **必需**): 输入的 `AnnData` 对象。`adata.raw` 必须包含原始计数数据。
    *   `pathway_dict` (dict, **必需**): `{pathway_name: gene_set}` 格式的通路基因集字典。
    *   `custom_genes` (list, 可选): 一个额外的自定义基因列表，将以 `'Custom'` 作为通路名称添加到 `pathway_dict` 中。
    *   `methods` (list of str, 默认值: `['AUCell', 'AddModule', 'ssGSEA', 'singscore']`): 要使用的打分方法列表。
        *   可选值: `'AUCell'`, `'AddModule'`, `'ssGSEA'`, `'singscore'`。
    *   `normalize` (bool, 默认值: `True`): 是否同时生成 Z-score 归一化的通路得分。如果为 `True`，则会为每个通路生成 `METHOD_PATHWAYNAME` (原始得分) 和 `METHOD_PATHWAYNAME_Z` (归一化得分) 两列。
    *   `min_genes` (int, 默认值: `2`): 通路中至少有多少个基因存在于 `adata.raw` 中才能进行打分。
*   **返回**:
    *   `adata`: 更新后的 `AnnData` 对象，其中 `adata.obs` 包含了通路得分。
    *   `method_stats`: 包含每种方法打分统计信息的字典。
*   **示例**:

    ```python
    import scGeneset as sci
    import scanpy as sc
    import os

    # 假设 adata 已经加载并包含 .raw 属性
    # adata = sci.load_data(adata_file='path/to/your/annotated_adata.h5ad')
    # 使用你的实际 adata 加载代码
    adata = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata.raw is None:
        adata.raw = adata.copy()

    pathway_dict = sci.load_pathway_genes('PCD')
    custom_genes = ['MMP10', 'CCL11', 'C4BPB', 'IL1RN', 'ANXA1']

    adata_scored, method_stats = sci.score_pathways(
        adata,
        pathway_dict,
        custom_genes=custom_genes,
        methods=['AUCell', 'AddModule'], # 仅使用 AUCell 和 AddModule 进行演示
        normalize=True
    )

    print("\nAUCell Custom 归一化得分 (前5个细胞):")
    print(adata_scored.obs['AUCell_Custom_Z'].head())
    ```

### `scGeneset.plot`

#### `plot_heatmap(data, method='AUCell', group=None, cell_type=None, z=True, pathway=None, custom_genes=None)`

*   **功能**: 绘制通路活性热图矩阵，展示通路得分在不同细胞类型和分组中的平均水平。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `method` (str, 默认值: `'AUCell'`): 要可视化的打分方法。
        *   可选值: `'AUCell'`, `'AddModule'`, `'ssGSEA'`, `'singscore'`。
    *   `group` (str or list of str, 可选): 要包含的分组名称。如果为 `None`，则包含所有分组。
    *   `cell_type` (str or list of str, 可选): 要包含的细胞类型名称。如果为 `None`，则包含所有细胞类型和“All Cells”汇总。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `pathway` (list of str, 可选): 要包含的通路名称列表。如果为 `None`，则包含所有可用通路。
    *   `custom_genes` (list, 可选): 一个自定义基因列表。如果 `pathway` 参数中包含 `'Custom'`，此参数用于识别 `Custom` 通路的基因。
*   **返回**:
    *   `fig`: `matplotlib.figure.Figure` 对象。
    *   `data_matrix`: 热图数据矩阵 `pandas.DataFrame`。
    *   `p_matrix`: p 值矩阵 `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    # 假设 adata_scored 已经通过 score_pathways 获得
    # adata_scored = ...
    # from your adata_scored (after running score_pathways)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()

    fig_heat, data_matrix, p_matrix = sci.plot_heatmap(
        adata_scored,
        method='AUCell',
        z=True,
        cell_type=None, # 所有细胞类型 + All Cells
        group=['HC', 'CD'] # 仅显示HC和CD组
    )
    fig_heat.savefig('heatmap_example.png', dpi=300, bbox_inches='tight')
    ```

#### `plot_violin(data, pathway, method='AUCell', group=None, cell_type=None, z=True, custom_genes=None)`

*   **功能**: 绘制指定通路的活性得分小提琴图，用于比较不同细胞类型或不同分组间的通路活性。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `pathway` (str, **必需**): 要绘制的通路名称。如果针对自定义基因集，传入 `'Custom'`。
    *   `method` (str, 默认值: `'AUCell'`): 要可视化的打分方法。
    *   `group` (str or list of str, 可选): 如果指定，则在 **这个分组** 内比较不同的细胞类型。可以是一个字符串或字符串列表。**不能与 `cell_type` 同时指定。**
    *   `cell_type` (str or list of str, 可选): 如果指定，则在 **这个细胞类型** 内比较不同的分组。可以是一个字符串或字符串列表。**不能与 `group` 同时指定。**
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `custom_genes` (list, 可选): 如果 `pathway` 是 `'Custom'`，此参数用于识别 `Custom` 通路的基因。
*   **返回**: `fig`: `matplotlib.figure.Figure` 对象。
*   **示例**:

    ```python
    import scGeneset as sci
    # adata_scored = ... (从 score_pathways 获得)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()
    custom_genes = ['MMP10', 'CCL11', 'C4BPB', 'IL1RN', 'ANXA1'] # 示例自定义基因列表

    # 示例1: 比较 'CD' 组内不同细胞类型中 'Custom' 通路的活性
    fig_violin1 = sci.plot_violin(
        adata_scored,
        pathway='Custom',
        method='AUCell',
        group='CD',         # 固定分组为 'CD'
        cell_type=None,     # 比较所有细胞类型
        z=True,
        custom_genes=custom_genes
    )
    fig_violin1.savefig('violin_custom_across_types_CD_example.png', dpi=300, bbox_inches='tight')

    # 示例2: 比较 'Paneth cells' 细胞类型内不同分组中 'Custom' 通路的活性
    fig_violin2 = sci.plot_violin(
        adata_scored,
        pathway='Custom',
        method='AUCell',
        group=None,         # 比较所有分组
        cell_type='Paneth cells', # 固定细胞类型为 'Paneth cells'
        z=True,
        custom_genes=custom_genes
    )
    fig_violin2.savefig('violin_custom_Paneth_cells_example.png', dpi=300, bbox_inches='tight')
    ```

#### `plot_umap(data, method='AUCell', z=True, top_n=5)`

*   **功能**: 绘制 UMAP 降维图，并在图周围的环形区域展示不同细胞类型中通路活性的平均水平和分组比例。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。`adata.obsm['X_umap']` 必须存在。
    *   `method` (str, 默认值: `'AUCell'`): 要可视化的打分方法。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `top_n` (int, 默认值: `5`): 在环形图的内圈中显示通路分数的前 `top_n` 个通路。
*   **返回**: `fig`: `matplotlib.figure.Figure` 对象。
*   **示例**:

    ```python
    import scGeneset as sci
    import scanpy as sc
    # adata_scored = ... (从 score_pathways 获得)
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()
    # 确保 UMAP 坐标已计算
    if 'X_umap' not in adata_scored.obsm:
        sc.pp.neighbors(adata_scored)
        sc.tl.umap(adata_scored)

    fig_umap = sci.plot_umap(
        adata_scored,
        method='AUCell',
        z=True,
        top_n=5
    )
    fig_umap.savefig('circular_umap_pathways_example.png', dpi=300, bbox_inches='tight')
    ```

#### `plot_corr_scatter(data, method='AUCell', group=None, cell_type=None, z=True, custom_genes=None)`

*   **功能**: 绘制自定义基因集得分与多个通路得分之间的散点图矩阵。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `group` (str, 可选): 要过滤的特定分组。如果为 `None`，则使用所有分组。
    *   `cell_type` (str, 可选): 要过滤的特定细胞类型。如果为 `None`，则使用所有细胞类型。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `custom_genes` (list, **必需**): 自定义基因列表。
*   **返回**: `fig`: `matplotlib.figure.Figure` 对象。
*   **示例**:

    ```python
    import scGeneset as sci
    import scanpy as sc
    # adata_scored = ... (从 score_pathways 获得)
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()
    custom_genes = ['MMP10', 'CCL11', 'C4BPB', 'IL1RN', 'ANXA1'] # 示例自定义基因列表

    fig_corr_scatter = sci.plot_corr_scatter(
        adata_scored,
        method='AUCell',
        group='CD', # 仅分析CD组
        cell_type='Paneth cells', # 仅分析Paneth cells
        z=True,
        custom_genes=custom_genes
    )
    fig_corr_scatter.savefig('custom_scatter_Paneth_CD_example.png', dpi=300, bbox_inches='tight')
    ```

#### `plot_corr_matrix(data, method='AUCell', group=None, cell_type=None, z=True, custom_genes=None)`

*   **功能**: 绘制通路间的相关性矩阵，并可视化自定义基因集与各个通路的相关性。该函数会调整共同基因集的影响。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。`adata.raw` 必须包含原始计数。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `group` (str, 可选): 要过滤的特定分组。如果为 `None`，则使用所有分组。
    *   `cell_type` (str, 可选): 要过滤的特定细胞类型。如果为 `None`，则使用所有细胞类型。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `custom_genes` (list, **必需**): 自定义基因列表。
*   **返回**:
    *   `fig`: `matplotlib.figure.Figure` 对象。
    *   `corr_matrix`: 相关性矩阵 `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    import scanpy as sc
    # adata_scored = ... (从 score_pathways 获得)
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()
    custom_genes = ['MMP10', 'CCL11', 'C4BPB', 'IL1RN', 'ANXA1'] # 示例自定义基因列表

    fig_corr_matrix, corr_matrix = sci.plot_corr_matrix(
        adata_scored,
        method='AUCell',
        group='CD', # 仅分析CD组
        cell_type='Plasma cells', # 仅分析Plasma cells
        z=True,
        custom_genes=custom_genes
    )
    fig_corr_matrix.savefig('correlation_custom_plasma_CD_example.png', dpi=300, bbox_inches='tight')
    corr_matrix.to_csv('correlation_matrix_plasma_CD_example.csv')
    ```

#### `plot_pca(data, method='AUCell', z=True, pathway=None)`

*   **功能**: 对通路得分执行 PCA (主成分分析)，并绘制 PCA 散点图和特征重要性条形图。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `pathway` (list of str, 可选): 要包含在 PCA 中的通路列表。如果为 `None`，则包含所有可用通路。
*   **返回**:
    *   `fig_scatter`: PCA 散点图 `matplotlib.figure.Figure` 对象。
    *   `fig_importance`: 特征重要性条形图 `matplotlib.figure.Figure` 对象。
    *   `loadings_df`: 载荷矩阵 `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    # adata_scored = ... (从 score_pathways 获得)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()

    fig_pca_scatter, fig_pca_importance, loadings_df = sci.plot_pca(
        adata_scored,
        method='AUCell',
        z=True,
        pathway=None # 使用所有通路进行PCA
    )
    fig_pca_scatter.savefig('pca_scatter_example.png', dpi=300, bbox_inches='tight')
    fig_pca_importance.savefig('pca_importance_example.png', dpi=300, bbox_inches='tight')
    loadings_df.to_csv('pca_loadings_example.csv')
    ```

### `scGeneset.analyze`

#### `calc_consistency(data, method='AUCell', group=None, z=True, custom_genes=None)`

*   **功能**: 计算自定义基因集得分与各个通路得分在不同细胞类型中的一致性（Spearman 相关性），并绘制点图。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `group` (str or list of str, 可选): 要过滤的特定分组。如果为 `None`，则使用所有分组。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `custom_genes` (list, **必需**): 自定义基因列表。
*   **返回**:
    *   `fig`: `matplotlib.figure.Figure` 对象。
    *   `consist_df`: 一致性结果 `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    # adata_scored = ... (从 score_pathways 获得)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()
    custom_genes = ['MMP10', 'CCL11', 'C4BPB', 'IL1RN', 'ANXA1'] # 示例自定义基因列表

    fig_consist, consist_df = sci.calc_consistency(
        adata_scored,
        method='AUCell',
        group=None, # 在所有分组中进行分析
        z=True,
        custom_genes=custom_genes
    )
    fig_consist.savefig('consistency_dotplot_example.png', dpi=300, bbox_inches='tight')
    consist_df.to_csv('consistency_results_example.csv', index=False)
    ```

#### `calc_wgcna(data, method='AUCell', group=None, cell_type=None, z=True, pathway=None, custom_genes=None)`

*   **功能**: 计算通路间的 WGCNA 共表达网络 (TOM 矩阵)，并绘制热图。该函数会调整共同基因集的影响。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。`adata.raw` 必须包含原始计数。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `group` (str, 可选): 要过滤的特定分组。如果为 `None`，则使用所有分组。
    *   `cell_type` (str, 可选): 要过滤的特定细胞类型。如果为 `None`，则使用所有细胞类型。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `pathway` (list of str, 可选): 要包含在 WGCNA 中的通路列表。如果为 `None`，则包含所有可用通路。
    *   `custom_genes` (list, 可选): 如果 `pathway` 参数中包含 `'Custom'`，此参数用于识别 `Custom` 通路的基因。
*   **返回**:
    *   `fig`: `matplotlib.figure.Figure` 对象。
    *   `tom_df`: TOM (拓扑重叠矩阵) `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    # adata_scored = ... (从 score_pathways 获得)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()

    fig_wgcna, tom_df = sci.calc_wgcna(
        adata_scored,
        method='AUCell',
        group=None, # 在所有分组中进行分析
        cell_type=None, # 在所有细胞类型中进行分析
        z=True,
        pathway=None # 使用所有通路
    )
    fig_wgcna.savefig('scwgcna_tom_heatmap_example.png', dpi=300, bbox_inches='tight')
    tom_df.to_csv('scwgcna_tom_matrix_example.csv')
    ```

#### `analyze_communication(data, method='AUCell', group=None, pathway=None, z=True)`

*   **功能**: 使用 LIANA (Ligand-Receptor Interaction Analysis) 进行细胞-细胞通讯分析，并生成热图和弦图。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法 (仅用于列名查找)。
    *   `group` (str, 可选): 要分析的特定分组。如果为 `None`，则分析所有细胞。
    *   `pathway` (str, 可选): 如果指定，则根据此通路的得分将细胞分为“高得分”和“低得分”两组进行分析。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分 (仅用于 `pathway` 参数)。
*   **返回**: 字典，包含每个分析组的结果 (例如 `'All'`, `'High'`, `'Low'`)，每个结果包含：
    *   `'liana_df'`: LIANA 结果 `pandas.DataFrame`。
    *   `'fig_heatmap'`: 热图 `matplotlib.figure.Figure` 对象。
    *   `'fig_chord'`: 弦图 `matplotlib.figure.Figure` 对象。
    *   `'interaction_matrix'`: 交互矩阵 `pandas.DataFrame`。
*   **注意**: 此功能需要安装 `liana` 库 (`pip install liana`)。
*   **示例**:

    ```python
    import scGeneset as sci
    # adata_scored = ... (从 score_pathways 获得)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()

    # 示例1: 所有细胞的细胞通讯
    comm_results_all = sci.analyze_communication(
        adata_scored,
        method='AUCell',
        group=None,
        pathway=None,
        z=True
    )
    for group_name, results in comm_results_all.items():
        results['fig_heatmap'].savefig(f'communication_heatmap_{group_name}_example.png', dpi=300, bbox_inches='tight')
        results['fig_chord'].savefig(f'communication_chord_{group_name}_example.png', dpi=300, bbox_inches='tight')

    # 示例2: 根据 'Intrinsic Apoptosis' 通路高低分组的细胞通讯
    comm_results_pathway_grouped = sci.analyze_communication(
        adata_scored,
        method='AUCell',
        group=None,
        pathway='Intrinsic Apoptosis', # 根据此通路得分高低分组
        z=True
    )
    for group_name, results in comm_results_pathway_grouped.items():
        results['fig_heatmap'].savefig(f'communication_heatmap_{group_name}_Intrinsic_Apoptosis_example.png', dpi=300, bbox_inches='tight')
        results['fig_chord'].savefig(f'communication_chord_{group_name}_Intrinsic_Apoptosis_example.png', dpi=300, bbox_inches='tight')
    ```

#### `analyze_pseudotime(data, method='AUCell', group=None, z=True, pathway=None, custom_genes=None, root_cell_type='Enterocytes', root_group_for_pseudotime=None)`

*   **功能**: 进行拟时序分析，计算细胞轨迹，并展示通路活性沿拟时序的趋势。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法 (用于通路得分列)。
    *   `group` (str, 可选): 要分析的特定分组。如果为 `None`，则分析所有细胞。**如果指定，`root_group_for_pseudotime` 将被忽略。**
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `pathway` (list of str, 可选): 要在热图和轨迹图中显示的通路列表。如果为 `None`，则包含所有可用通路。
    *   `custom_genes` (dict, 可选): 自定义基因集字典。如果 `pathway` 参数中包含 `'Custom'`，此参数用于识别 `Custom` 通路的基因。
    *   `root_cell_type` (str, 默认值: `'Enterocytes'`): 用于识别拟时序根细胞的细胞类型。
    *   `root_group_for_pseudotime` (str, 可选): 如果 `group` 参数为 `None`，此参数指定在哪个分组中寻找 `root_cell_type` 来确定根细胞。如果此参数和 `group` 都为 `None`，则使用 `data.obs['group']` 的第一个类别。
*   **返回**: 字典，包含：
    *   `'fig_umap'`: 拟时序 UMAP `matplotlib.figure.Figure` 对象。
    *   `'fig_heatmap'`: 拟时序通路热图 `matplotlib.figure.Figure` 对象。
    *   `'fig_trajectory'`: 拟时序轨迹图 `matplotlib.figure.Figure` 对象。
    *   `'pseudotime'`: 拟时序值 `pandas.Series`。
*   **示例**:

    ```python
    import scGeneset as sci
    import scanpy as sc
    # adata_scored = ... (从 score_pathways 获得)
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()
    # 确保 UMAP 坐标已计算
    if 'X_umap' not in adata_scored.obsm:
        sc.pp.neighbors(adata_scored)
        sc.tl.umap(adata_scored)

    pseudo_results = sci.analyze_pseudotime(
        adata_scored,
        method='AUCell',
        group=None, # 在所有细胞中进行拟时序
        z=True,
        pathway=['Intrinsic Apoptosis', 'MAPK signaling pathway'], # 指定要显示的通路
        root_cell_type='Enterocytes', # 假设 'Enterocytes' 是你的起始细胞类型
        root_group_for_pseudotime='HC' # 假设在'HC'组中寻找根细胞
    )
    pseudo_results['fig_umap'].savefig('pseudotime_umap_example.png', dpi=300, bbox_inches='tight')
    pseudo_results['fig_heatmap'].savefig('pseudotime_heatmap_example.png', dpi=300, bbox_inches='tight')
    if pseudo_results['fig_trajectory'] is not None:
        pseudo_results['fig_trajectory'].savefig('pseudotime_trajectory_example.png', dpi=300, bbox_inches='tight')
    pseudo_results['pseudotime'].to_csv('pseudotime_values_example.csv')
    ```

#### `recluster_cells(data, pathway, method='AUCell', group=None, cell_type='Plasma cells', z=True, cluster='high')`

*   **功能**: 根据特定通路的得分（高或低），对指定细胞类型进行子聚类，并绘制 UMAP 和 marker 基因气泡图。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `pathway` (str, **必需**): 要根据其得分进行筛选的通路名称。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `group` (str, 可选): 要过滤的特定分组。如果为 `None`，则使用所有分组。
    *   `cell_type` (str, 默认值: `'Plasma cells'`): 要进行子聚类的细胞类型。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
    *   `cluster` (str, 默认值: `'high'`): 选择高得分细胞 (`'high'`) 或低得分细胞 (`'low'`) 进行子聚类。
        *   可选值: `'high'`, `'low'`。
*   **返回**: 字典，包含：
    *   `'fig_umap'`: 子聚类 UMAP `matplotlib.figure.Figure` 对象。
    *   `'fig_bubble'`: marker 基因气泡图 `matplotlib.figure.Figure` 对象。
    *   `'adata_sub'`: 子聚类后的 `AnnData` 对象。
    *   `'hvgs'`: 高度可变基因 `pandas.DataFrame`。
    *   `'markers'`: marker 基因 `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    import scanpy as sc
    # adata_scored = ... (从 score_pathways 获得)
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()

    reclust_results = sci.recluster_cells(
        adata_scored,
        pathway='Intrinsic Apoptosis',
        method='AUCell',
        group=None,
        cell_type='Plasma cells',
        z=True,
        cluster='high'
    )
    reclust_results['fig_umap'].savefig('plasma_subcluster_umap_example.png', dpi=300, bbox_inches='tight')
    if reclust_results['fig_bubble'] is not None:
        reclust_results['fig_bubble'].savefig('plasma_subcluster_bubble_example.png', dpi=300, bbox_inches='tight')
    reclust_results['hvgs'].to_csv('plasma_subcluster_hvgs_example.csv', index=False)
    reclust_results['markers'].to_csv('plasma_subcluster_markers_example.csv', index=False)
    ```

#### `analyze_differential(data, pathway, method='AUCell', group=None, cell_type='Plasma cells', z=True)`

*   **功能**: 对特定通路高得分细胞与低得分细胞之间进行差异基因分析，并绘制火山图。
*   **参数**:
    *   `data` (`AnnData`, **必需**): 输入的 `AnnData` 对象。
    *   `pathway` (str, **必需**): 要根据其得分进行分组的通路名称。
    *   `method` (str, 默认值: `'AUCell'`): 要使用的打分方法。
    *   `group` (str, 可选): 要过滤的特定分组。如果为 `None`，则使用所有分组。
    *   `cell_type` (str, 默认值: `'Plasma cells'`): 要进行差异分析的细胞类型。
    *   `z` (bool, 默认值: `True`): 是否使用 Z-score 归一化的通路得分。
*   **返回**: 字典，包含：
    *   `'fig_volcano'`: 火山图 `matplotlib.figure.Figure` 对象。
    *   `'deg_df'`: 差异表达基因 `pandas.DataFrame`。
*   **示例**:

    ```python
    import scGeneset as sci
    # adata_scored = ... (从 score_pathways 获得)
    import scanpy as sc
    adata_scored = sc.read_h5ad('./单细胞输出/final_data_with_scores.h5ad')
    if adata_scored.raw is None:
        adata_scored.raw = adata_scored.copy()

    diff_results = sci.analyze_differential(
        adata_scored,
        pathway='Custom',
        method='AUCell',
        group=None,
        cell_type='Plasma cells',
        z=True
    )
    diff_results['fig_volcano'].savefig('volcano_plasma_cells_custom_example.png', dpi=300, bbox_inches='tight')
    diff_results['deg_df'].to_csv('dea_plasma_cells_custom_example.csv', index=False)
    ```

### `scGeneset.utils`

此模块主要包含内部配置和辅助函数，通常不需要直接调用。
*   `setup_plotting()`: 设置全局 `matplotlib` 绘图参数。在库导入时自动调用。
*   `simplify_name(pathway, method='')`: 简化通路名称以便显示。
*   `format_pval(p)`: 格式化 p 值字符串。
*   `calc_gene_score_subset(...)`: 在细胞子集中计算基因集得分（内部辅助函数）。
*   `get_residuals(...)`: 计算 OLS 回归残差（内部辅助函数）。
*   `GLOBAL_PALETTE`, `GROUP_COLORS`, `CMAP_POS`, `CMAP_NEG_POS`, `CMAP_GREEN_RED`, `LEGEND_KWARGS`, `N_JOBS`, `PATHWAY_FILE`: 全局常量和配色方案。

---

## 完整运行示例


```python
# ==============================================================================
# scGeneset - 完整分析流程示例
# ==============================================================================

import sys
import os
import scanpy as sc
import anndata
import matplotlib.pyplot as plt # 导入 matplotlib 用于关闭所有图表
import numpy as np 

# ===== 添加当前目录到Python路径 (如果作为独立脚本运行，而非安装包) =====
# current_dir = os.path.dirname(os.path.abspath(__file__))
# sys.path.insert(0, current_dir)
# 在 Jupyter/Colab 环境中，可以直接导入，无需此步
# import scGeneset as sci # 如果已经安装为 scGeneset
import scGeneset as sci 
print(f"库位置: {sci.__file__}")
print(f"当前工作目录: {os.getcwd()}")

# ===========================
# 0. 准备输出目录
# ===========================
base_output_dir = './单细胞输出_scGeneset_example/'
os.makedirs(base_output_dir, exist_ok=True)

# ===========================
# 1. 数据加载与通路打分
# ===========================
# 定义一个 AnnData 对象的路径作为示例输入。
# 你需要替换为你的实际文件路径，或者使用 load_data 的文件列表模式。
# 这里假设你已经有一个名为 'annotated_adata.h5ad' 的文件
# 该文件应包含 .X (原始表达量), .obs (有'final_annotation'和'group'列), .var
# 并且最好有 .raw 属性存储原始计数
example_adata_path = os.path.join(base_output_dir, 'annotated_adata.h5ad')

# --- 假数据生成 (如果你的 annotated_adata.h5ad 不存在，运行此段生成) ---
if not os.path.exists(example_adata_path):
    print("生成示例 AnnData 文件...")
    adata_temp = sc.AnnData(X=sc.external.pp.downsample_counts(
        np.random.rand(200, 1000) * 100, 1000, copy=True
    ).X, dtype=np.float32)
    adata_temp.obs_names = [f'cell_{i}' for i in range(adata_temp.n_obs)]
    adata_temp.var_names = [f'gene_{i}' for i in range(adata_temp.n_vars)]
    adata_temp.obs['final_annotation'] = np.random.choice(
        ['T_cell', 'B_cell', 'Plasma cells', 'Paneth cells', 'Enterocytes'], adata_temp.n_obs)
    adata_temp.obs['group'] = np.random.choice(['HC', 'UC', 'CD'], adata_temp.n_obs)
    # 模拟原始数据
    adata_temp.raw = adata_temp.copy()
    adata_temp.write(example_adata_path)
    print(f"示例 AnnData 文件已保存到: {example_adata_path}")
# --- 假数据生成结束 ---

adata = sci.load_data(adata_file=example_adata_path)

# 确保 adata.raw 存在且是 AnnData 对象
if adata.raw is None:
    print("⚠️  警告: adata.raw 为 None，将使用当前数据作为 raw")
    adata.raw = adata.copy()
elif not isinstance(adata.raw, anndata.AnnData):
    print(f"⚠️  警告: adata.raw 类型为 {type(adata.raw)}，正在转换为 AnnData...")
    if hasattr(adata.raw, 'to_adata'):
        adata.raw = adata.raw.to_adata()
        print("✓ 成功转换 adata.raw 为 AnnData 对象")
    else:
        print("⚠️  无法转换，将使用当前数据作为 raw")
        adata.raw = adata.copy()
else:
    print("✓ adata.raw 已经是 AnnData 对象")


pathway_dict = sci.load_pathway_genes(pathway_file="PCD")

CUSTOM_GENES = ['MMP10', 'CCL11', 'C4BPB', 'IL1RN', 'ANXA1', 'PI3',
                'VDR', 'MMP3', 'TIMP1', 'S100A8', 'S100A9']
pathway_dict['Custom'] = set(CUSTOM_GENES)

# 为了演示细胞通讯，LIANA 需要更多基因集，这里添加一些假基因集
if len(pathway_dict) < 20: # 确保有足够多的通路进行WGCNA和通讯分析
    for i in range(20 - len(pathway_dict)):
        pathway_dict[f'FakePathway_{i}'] = set(adata.var_names[np.random.choice(adata.n_vars, 20, replace=False)].tolist())

print(f"✓ 通路加载完成，共 {len(pathway_dict)} 个通路")
print(f"  包含预定义通路 + Custom通路")

print("\n" + "=" * 60)
print("步骤3: 计算通路得分 (这可能需要几分钟...)")
print("=" * 60)

adata, method_stats = sci.score_pathways(
    adata,
    pathway_dict,
    methods=['AUCell'], # 仅使用 AUCell 演示，避免运行时间过长
    normalize=True,
)

print(f"\n✓ 打分完成！")
for method, stats in method_stats.items():
    score_cols = [col for col in adata.obs.columns if col.startswith(f'{method}_')]
    print(f"  - {method}: {len(score_cols)} 个通路得分列")

output_scored_path = os.path.join(base_output_dir, 'final_data_with_scores.h5ad')
adata.write(output_scored_path)
print(f"✓ 打分数据已保存到: {output_scored_path}")

print("\n" + "=" * 60)
print("步骤4: 快速验证打分结果")
print("=" * 60)
col_name = 'AUCell_Custom'
if col_name in adata.obs.columns:
    print(f"\n{col_name} (前5个细胞):")
    print(adata.obs[col_name].head())
print("\n✅ 打分步骤完成！")

# ===========================
# 2. 分析与可视化
# ===========================

print("\n" + "=" * 70)
print("步骤5: 生成 AUCell 分析图表")
print("=" * 70)

method = 'AUCell'
analysis_output_dir = os.path.join(base_output_dir, method)
os.makedirs(analysis_output_dir, exist_ok=True)


# --- 4.1 热图 ---
print(f"\n[{method}] 生成热图...")
fig_heat, data_matrix, p_matrix = sci.plot_heatmap(
    adata,
    method=method,
    z=True,
    cell_type=None, # 所有细胞类型 + All Cells
    group=None  # 所有分组
)
fig_heat.savefig(os.path.join(analysis_output_dir, 'pathway_heatmap_matrix.png'),
                 dpi=300, bbox_inches='tight')
data_matrix.to_csv(os.path.join(analysis_output_dir, 'heatmap_data.csv'))
p_matrix.to_csv(os.path.join(analysis_output_dir, 'heatmap_pvalues.csv'))
plt.close(fig_heat) # 关闭图表释放内存
print(f"  ✓ 热图已保存")

# --- 4.2 小提琴图 ---
print(f"[{method}] 生成小提琴图...")
fig_violin1 = sci.plot_violin(
    adata,
    pathway='Custom',
    method=method,
    group='CD', # 假设 'CD' 存在
    cell_type=None,
    z=True,
    custom_genes=CUSTOM_GENES
)
fig_violin1.savefig(os.path.join(analysis_output_dir, 'violin_custom_across_types_CD.png'),
                    dpi=300, bbox_inches='tight')
plt.close(fig_violin1)

fig_violin2 = sci.plot_violin(
    adata,
    pathway='Custom',
    method=method,
    group=None,
    cell_type='Paneth cells', # 假设 'Paneth cells' 存在
    z=True,
    custom_genes=CUSTOM_GENES
)
fig_violin2.savefig(os.path.join(analysis_output_dir, 'violin_custom_Paneth_cells.png'),
                    dpi=300, bbox_inches='tight')
plt.close(fig_violin2)
print(f"  ✓ 小提琴图已保存")

# --- 4.3 环形UMAP ---
print(f"[{method}] 生成环形UMAP...")
# 确保 UMAP 坐标已计算
if 'X_umap' not in adata.obsm:
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

fig_umap = sci.plot_umap(
    adata,
    method=method,
    z=True,
    top_n=5
)
fig_umap.savefig(os.path.join(analysis_output_dir, 'circular_umap_pathways.png'),
                 dpi=300, bbox_inches='tight')
plt.close(fig_umap)
print(f"  ✓ UMAP已保存")

# --- 4.4 相关性散点图 ---
print(f"[{method}] 生成Custom相关性散点图...")
for group_val in ['HC', 'CD']: # 假设 'HC', 'CD' 存在
    fig_corr_scatter = sci.plot_corr_scatter(
        adata,
        method=method,
        group=group_val,
        cell_type='Paneth cells', # 假设 'Paneth cells' 存在
        z=True,
        custom_genes=CUSTOM_GENES
    )
    fig_corr_scatter.savefig(
        os.path.join(analysis_output_dir, f'Custom_scatter_Paneth_{group_val}.png'),
        dpi=300, bbox_inches='tight'
    )
    plt.close(fig_corr_scatter)
print(f"  ✓ 相关性散点图已保存")

# --- 4.5 相关性矩阵 ---
print(f"[{method}] 生成Custom相关性矩阵...")
for group_val in ['HC', 'UC', 'CD']: # 假设 'HC', 'UC', 'CD' 存在
    fig_corr_matrix, corr_matrix = sci.plot_corr_matrix(
        adata,
        method=method,
        group=group_val,
        cell_type='Plasma cells', # 假设 'Plasma cells' 存在
        z=True,
        custom_genes=CUSTOM_GENES
    )
    fig_corr_matrix.savefig(
        os.path.join(analysis_output_dir, f'correlation_Custom_Plasma_{group_val}.png'),
        dpi=300, bbox_inches='tight'
    )
    corr_matrix.to_csv(
        os.path.join(analysis_output_dir, f'correlation_matrix_Plasma_{group_val}.csv')
    )
    plt.close(fig_corr_matrix)
print(f"  ✓ 相关性矩阵已保存")

# --- 4.6 PCA分析 ---
print(f"[{method}] 生成PCA图...")
fig_pca_scatter, fig_pca_importance, loadings_df = sci.plot_pca(
    adata,
    method=method,
    z=True,
    pathway=None
)
fig_pca_scatter.savefig(os.path.join(analysis_output_dir, 'pca_scatter.png'),
                        dpi=300, bbox_inches='tight')
fig_pca_importance.savefig(os.path.join(analysis_output_dir, 'pca_importance.png'),
                           dpi=300, bbox_inches='tight')
loadings_df.to_csv(os.path.join(analysis_output_dir, 'pca_loadings.csv'))
plt.close(fig_pca_scatter)
plt.close(fig_pca_importance)
print(f"  ✓ PCA分析已保存")

# --- 4.7 一致性分析 ---
print(f"[{method}] 生成一致性点图...")
fig_consist, consist_df = sci.calc_consistency(
    adata,
    method=method,
    group=None,
    z=True,
    custom_genes=CUSTOM_GENES
)
fig_consist.savefig(os.path.join(analysis_output_dir, 'consistency_dotplot.png'),
                    dpi=300, bbox_inches='tight')
consist_df.to_csv(os.path.join(analysis_output_dir, 'consistency_results.csv'),
                  index=False)
plt.close(fig_consist)
print(f"  ✓ 一致性分析已保存")

# --- 4.8 WGCNA分析 ---
print(f"[{method}] 运行WGCNA分析...")
fig_wgcna, tom_df = sci.calc_wgcna(
    adata,
    method=method,
    group=None,
    cell_type=None,
    z=True,
    pathway=None
)
fig_wgcna.savefig(os.path.join(analysis_output_dir, 'scwgcna_tom_heatmap.png'),
                  dpi=300, bbox_inches='tight')
tom_df.to_csv(os.path.join(analysis_output_dir, 'scwgcna_tom_matrix.csv'))
plt.close(fig_wgcna)
print(f"  ✓ WGCNA分析已保存")

# --- 4.9 细胞通讯分析 ---
print(f"[{method}] 运行细胞通讯分析...")
try:
    import liana as li
    # 所有细胞的通讯
    comm_results_all = sci.analyze_communication(
        adata,
        method=method,
        group=None,
        pathway=None,
        z=True
    )
    for group_name, results in comm_results_all.items():
        results['liana_df'].to_csv(
            os.path.join(analysis_output_dir, f'liana_results_{group_name}.csv'),
            index=False
        )
        results['fig_heatmap'].savefig(
            os.path.join(analysis_output_dir, f'communication_heatmap_{group_name}.png'),
            dpi=300, bbox_inches='tight'
        )
        results['fig_chord'].savefig(
            os.path.join(analysis_output_dir, f'communication_chord_{group_name}.png'),
            dpi=300, bbox_inches='tight'
        )
        results['interaction_matrix'].to_csv(
            os.path.join(analysis_output_dir, f'interaction_matrix_{group_name}.csv')
        )
        plt.close(results['fig_heatmap'])
        plt.close(results['fig_chord'])

    # 按 'Intrinsic Apoptosis' 通路分组
    comm_results_pathway_grouped = sci.analyze_communication(
        adata,
        method=method,
        group=None,
        pathway='Intrinsic Apoptosis', # 假设 'Intrinsic Apoptosis' 通路存在
        z=True
    )
    for group_name, results in comm_results_pathway_grouped.items():
        results['liana_df'].to_csv(
            os.path.join(analysis_output_dir, f'liana_results_{group_name}_Intrinsic_Apoptosis.csv'),
            index=False
        )
        results['fig_heatmap'].savefig(
            os.path.join(analysis_output_dir, f'communication_heatmap_{group_name}_Intrinsic_Apoptosis.png'),
            dpi=300, bbox_inches='tight'
        )
        results['fig_chord'].savefig(
            os.path.join(analysis_output_dir, f'communication_chord_{group_name}_Intrinsic_Apoptosis.png'),
            dpi=300, bbox_inches='tight'
        )
        plt.close(results['fig_heatmap'])
        plt.close(results['fig_chord'])

    print(f"  ✓ 细胞通讯分析已保存")

except ImportError:
    print("  ⚠️ LIANA 未安装，跳过细胞通讯分析。请运行 `pip install liana` 安装。")


# --- 4.10 拟时序分析 ---
print(f"[{method}] 运行拟时序分析...")
# 确保 diffmap 已计算
sc.tl.diffmap(adata)
pseudo_results = sci.analyze_pseudotime(
    adata,
    method=method,
    group=None,
    z=True,
    pathway=None,
    root_cell_type='Enterocytes', # 假设 'Enterocytes' 存在
    root_group_for_pseudotime='HC' # 假设 'HC' 存在
)

pseudo_results['fig_umap'].savefig(
    os.path.join(analysis_output_dir, 'pseudotime_umap.png'),
    dpi=300, bbox_inches='tight'
)
pseudo_results['fig_heatmap'].savefig(
    os.path.join(analysis_output_dir, 'pseudotime_heatmap.png'),
    dpi=300, bbox_inches='tight'
)
if pseudo_results['fig_trajectory'] is not None:
    pseudo_results['fig_trajectory'].savefig(
        os.path.join(analysis_output_dir, 'pseudotime_trajectory.png'),
        dpi=300, bbox_inches='tight'
    )
pseudo_results['pseudotime'].to_csv(
    os.path.join(analysis_output_dir, 'pseudotime_values.csv')
)
plt.close(pseudo_results['fig_umap'])
plt.close(pseudo_results['fig_heatmap'])
if pseudo_results['fig_trajectory'] is not None:
    plt.close(pseudo_results['fig_trajectory'])
print(f"  ✓ 拟时序分析已保存")

# --- 4.11 Plasma细胞重聚类 ---
print(f"[{method}] Plasma细胞重聚类...")
reclust_results = sci.recluster_cells(
    adata,
    pathway='Intrinsic Apoptosis', # 假设 'Intrinsic Apoptosis' 通路存在
    method=method,
    group=None,
    cell_type='Plasma cells', # 假设 'Plasma cells' 存在
    z=True,
    cluster='high'
)

reclust_results['fig_umap'].savefig(
    os.path.join(analysis_output_dir, 'plasma_subcluster_umap.png'),
    dpi=300, bbox_inches='tight'
)
if reclust_results['fig_bubble'] is not None:
    reclust_results['fig_bubble'].savefig(
        os.path.join(analysis_output_dir, 'plasma_subcluster_bubble.png'),
        dpi=300, bbox_inches='tight'
    )
reclust_results['hvgs'].to_csv(
    os.path.join(analysis_output_dir, 'plasma_subcluster_hvgs.csv'),
    index=False
)
reclust_results['markers'].to_csv(
    os.path.join(analysis_output_dir, 'plasma_subcluster_markers.csv'),
    index=False
)
plt.close(reclust_results['fig_umap'])
if reclust_results['fig_bubble'] is not None:
    plt.close(reclust_results['fig_bubble'])
print(f"  ✓ Plasma细胞重聚类已保存")

# --- 4.12 差异分析 ---
print(f"[{method}] 差异表达分析...")
diff_results = sci.analyze_differential(
    adata,
    pathway='Custom',
    method=method,
    group=None,
    cell_type='Plasma cells', # 假设 'Plasma cells' 存在
    z=True
)

diff_results['fig_volcano'].savefig(
    os.path.join(analysis_output_dir, 'volcano_Plasma_cells_Custom.png'),
    dpi=300, bbox_inches='tight'
)
diff_results['deg_df'].to_csv(
    os.path.join(analysis_output_dir, 'dea_Plasma_cells_Custom.csv'),
    index=False
)
plt.close(diff_results['fig_volcano'])
print(f"  ✓ 差异分析已保存")

print(f"\n✅ {method} 分析完成!")

print("\n" + "=" * 70)
print("✅ 所有分析步骤完成！")
print("=" * 70)
print(f"\n结果保存在: {base_output_dir}{method}/")
