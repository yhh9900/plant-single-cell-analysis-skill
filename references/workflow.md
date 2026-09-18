# FASTQ 到 Harmony 注释流程

## 1. FASTQ 与样本信息检查

建立样本表，每一行对应一个文库，至少包括：样本名、FASTQ 路径、生物学分组、技术批次、建库平台、`scRNA` 或 `snRNA`、参考基因组版本和预期读段结构。

先确认哪条 read 保存细胞条形码和 UMI、哪条 read 保存转录本序列。短的 barcode/UMI read 不是无用 read，不能因为长度短而删除。可用 FastQC 和 MultiQC 汇总 read 数、长度、碱基质量、接头和重复率。

## 2. 从 FASTQ 生成表达矩阵

根据经过确认的建库化学选择计数工具：

- 10x Chromium：使用兼容版本的 Cell Ranger，或参数经过核实的 STARsolo。
- 其他液滴或组合索引平台：优先采用厂商流程，或使用已经验证的 STARsolo、kallisto|bustools 参数。
- 单核 RNA 测序：在平台支持时计入内含子 reads，并记录该设置。

参考索引使用选定版本的 genome FASTA 与 GTF/GFF 构建。保存命令、软件版本、参考文件版本、校验值、chemistry 参数和是否计入内含子。每个样本先独立生成矩阵并检查，之后再合并。

## 3. 单样本质量控制

每个样本分别创建 Seurat 对象，计算检测基因数、UMI 数、线粒体、叶绿体和核糖体比例。线粒体与叶绿体比例必须分开保存和筛选。基因集合或匹配规则应从当前物种和当前版本的 GTF/GFF 得到，不能直接照搬拟南芥正则表达式到其他植物。

根据每个样本的分布和生物学背景设置阈值，不机械复制统一阈值。视平台和可用软件评估空液滴、环境 RNA 和双细胞。输出每个样本每一步过滤前后的细胞数。

## 4. 标准化与 Harmony

保留原始 counts 不变。全部样本采用一致的标准化分支：

- `LogNormalize`、高变基因、`ScaleData`、PCA；或
- SCTransform 及与之匹配的 Harmony 流程。

Harmony 的变量必须代表技术批次或样本来源。若技术批次与生物学处理完全混杂，应报告混杂，不能直接把处理组当作批次消除。

PCA 后运行 Harmony，并基于 Harmony embedding 运行邻居图、聚类和 UMAP：

```r
obj <- RunHarmony(obj, group.by.vars = "batch", reduction = "pca")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_harmony_dims)
obj <- FindClusters(obj, resolution = cluster_resolution)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_harmony_dims)
```

维度数和 resolution 根据拐点图、聚类稳定性和生物学可解释性确定，不固定照抄某个项目的数值。保留未校正 PCA 与 Harmony 两套 reduction，比较批次混合和真实生物信号是否同时合理。

## 5. 先算 marker，再做注释

将当前聚类结果设为 identity，为每个 cluster 计算正向 marker。按校正后显著性、效应量和表达比例差异排序，导出：

- `markers_all.csv`：全部 marker 及 cluster、效应量、校正 P 值、表达比例；
- `markers_top10_by_cluster.csv`：每个 cluster 排名前 10 的 marker；
- Top10 marker 的 DotPlot 或 heatmap。

必须先完成这些结果，再填写细胞类型注释。

## 6. 匹配拟南芥基因名

表达矩阵、Seurat 对象和差异表达统计始终保留原物种基因 ID。拟南芥匹配只作为 marker 解释层加入，不直接替换矩阵行名，也不把多个同源基因的表达量静默合并。

- 拟南芥数据：把当前注释中的 feature ID 对应到 TAIR locus ID 和标准 symbol，同时保留两列。
- 其他植物：使用当前项目的同源基因表，把原物种 marker 对应到拟南芥 locus ID 和标准 symbol。
- 一对一同源关系可直接用于 marker 解释。
- 一对多、多对一和多对多关系全部保留，并在 `orthology_type` 中标记；棉花等多倍体的 homeolog 不能只保留一个。

使用 `scripts/map_markers_to_arabidopsis.R` 生成带拟南芥名称的完整 marker 表和 Top10 marker 表。注释时同时查看原物种 marker、拟南芥名称、同源关系类型和表达特异性。

## 7. 注释依据与命名

对每个 cluster 记录：基础注释、原物种支持 marker、拟南芥对应名称、冲突 marker、同源关系类型、可信度和参考来源。跨物种注释若依赖同源基因，应在结果中明确说明属于同源推断。

基础注释全部完成后执行重复名称规则：只有重复的细胞类型前面加一个区分该 cluster 的 marker，格式为 `Marker+ cell type`；不重复的细胞类型只保留细胞类型名称。

## 8. 可重复性

分别保存导入矩阵、QC 后、Harmony 后、聚类后和最终注释后的对象。记录随机种子、软件版本、参数、参考版本和人工注释依据。下游分析不覆盖原始 FASTQ 和原始表达矩阵。
