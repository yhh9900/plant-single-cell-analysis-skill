# 第一阶段：从每个样本的表达矩阵开始，完成 QC、Harmony、聚类、UMAP 和 Top10 marker。
# FASTQ 到表达矩阵的计数命令由建库平台决定，应先按照 references/workflow.md 完成。
#
# 运行方式：
# Rscript 01_harmony_cluster_markers.R 项目目录 sample_metadata.csv 线粒体基因正则表达式 叶绿体基因正则表达式
#
# sample_metadata.csv 每行对应一个样本，固定包含：
# sample,matrix_dir,batch,group,min_features,max_features,max_counts,max_percent_mt,max_percent_cp
# matrix_dir 指向 Cell Ranger、STARsolo 或其他计数工具输出的 10x 格式矩阵目录。
# 线粒体和叶绿体阈值应根据每个样本的分布提前填写，不在脚本中自动猜测。

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# 固定随机种子，使 PCA、UMAP 等含随机过程的结果可以重复。
set.seed(2026)

# 从命令行读取项目目录、样本表，以及当前植物物种的线粒体和叶绿体基因匹配规则。
args <- commandArgs(trailingOnly = TRUE)
project_dir <- args[1]
sample_metadata_file <- args[2]
mitochondrial_pattern <- args[3]
chloroplast_pattern <- args[4]

# 集中设置需要根据项目调整的分析参数。
# n_pcs 是 PCA 计算的维度数；n_harmony_dims 是建邻居图和 UMAP 使用的 Harmony 维度数。
# cluster_resolution 控制聚类粒度，需要结合稳定性和 marker 结果判断。
n_pcs <- 50
n_harmony_dims <- 30
cluster_resolution <- 0.6

# 建立输出目录，并读取已经准备好的样本信息表。
output_dir <- file.path(project_dir, "harmony_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sample_metadata <- read.csv(sample_metadata_file, stringsAsFactors = FALSE)

# 每个样本独立读取矩阵、建立 Seurat 对象并执行 QC。
# 独立 QC 可以避免某个样本的测序深度分布掩盖另一个样本的问题。
object_list <- lapply(seq_len(nrow(sample_metadata)), function(i) {
  sample_info <- sample_metadata[i, ]

  # 读取 gene-by-cell 稀疏矩阵，并使用样本名作为 project 标识。
  counts <- Read10X(data.dir = sample_info$matrix_dir)
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_info$sample,
    min.cells = 3
  )

  # 给细胞条形码添加样本前缀，避免合并不同样本时出现重名。
  obj <- RenameCells(obj, add.cell.id = sample_info$sample)

  # 写入 Harmony 和后续分组比较需要的样本元数据。
  obj$sample <- sample_info$sample
  obj$batch <- sample_info$batch
  obj$group <- sample_info$group

  # 线粒体和叶绿体 reads 代表不同的植物细胞器来源，因此分别计算比例。
  # 正则表达式应根据当前参考注释的 feature 名称提前确定。
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = mitochondrial_pattern
  )
  obj[["percent.cp"]] <- PercentageFeatureSet(
    obj,
    pattern = chloroplast_pattern
  )

  # 使用样本表中已经确定的阈值过滤低质量细胞和异常高计数细胞。
  obj <- subset(
    obj,
    subset = nFeature_RNA >= sample_info$min_features &
      nFeature_RNA <= sample_info$max_features &
      nCount_RNA <= sample_info$max_counts &
      percent.mt <= sample_info$max_percent_mt &
      percent.cp <= sample_info$max_percent_cp
  )

  obj
})

# 合并所有样本。细胞名称已经带样本前缀，因此合并后仍可追溯到原始文库。
combined <- Reduce(function(x, y) merge(x, y), object_list)

# 使用 LogNormalize 分支完成标准化、高变基因选择、缩放和 PCA。
# 原始 counts 保留在 RNA assay 中，marker 计算仍可追溯到原始计数。
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 3000)
combined <- ScaleData(combined, features = VariableFeatures(combined))
combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = n_pcs
)

# 输出 PCA 拐点图，便于检查设置的 Harmony 维度是否合理。
pca_elbow <- ElbowPlot(combined, ndims = n_pcs)
ggsave(
  file.path(output_dir, "pca_elbow.png"),
  pca_elbow,
  width = 6,
  height = 4,
  dpi = 300
)

# 使用技术批次 batch 运行 Harmony。
# Harmony 只校正低维表示，不修改原始表达矩阵。
combined <- RunHarmony(
  combined,
  group.by.vars = "batch",
  reduction = "pca"
)

# 邻居图、聚类和 UMAP 全部基于 Harmony embedding，保证流程口径一致。
combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:n_harmony_dims
)
combined <- FindClusters(
  combined,
  resolution = cluster_resolution
)
combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:n_harmony_dims
)

# 输出按 cluster 和 sample 着色的检查图，用于同时查看聚类结构和批次混合情况。
cluster_umap <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
sample_umap <- DimPlot(combined, reduction = "umap", group.by = "sample")
ggsave(file.path(output_dir, "umap_clusters_check.png"), cluster_umap, width = 7, height = 6, dpi = 300)
ggsave(file.path(output_dir, "umap_samples_check.png"), sample_umap, width = 7, height = 6, dpi = 300)

# 在任何细胞类型注释之前计算所有 cluster marker。
Idents(combined) <- "seurat_clusters"
markers_all <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
)

# 先按 cluster 分组，再按平均 log2 fold change 从高到低排列，提取每群前 10 个 marker。
markers_top10 <- markers_all %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj) %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

# 保存完整 marker 表和每群 Top10 marker 表。
# 下一步先运行 map_markers_to_arabidopsis.R 增加拟南芥 locus ID 和标准基因名，
# 然后再结合原物种 marker 与拟南芥同源 marker 进行人工注释。
write.csv(
  markers_all,
  file.path(output_dir, "markers_all.csv"),
  row.names = FALSE
)
write.csv(
  markers_top10,
  file.path(output_dir, "markers_top10_by_cluster.csv"),
  row.names = FALSE
)

# 用 Top10 marker 绘制 DotPlot，快速检查每个 marker 在所有 cluster 中的特异性。
top10_genes <- unique(markers_top10$gene)
top10_dotplot <- DotPlot(combined, features = top10_genes) +
  RotatedAxis()
ggsave(
  file.path(output_dir, "markers_top10_dotplot.png"),
  top10_dotplot,
  width = 14,
  height = 9,
  dpi = 300
)

# 保存注释前对象。下一阶段读取该对象并应用人工确认的注释表。
saveRDS(
  combined,
  file.path(output_dir, "seurat_harmony_clustered_before_annotation.rds")
)

# 保存软件版本，保证分析可以复现。
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
