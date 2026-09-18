# 第二阶段：读取 Top10 marker 审核后填写的注释表，生成最终注释和编号式 UMAP。
#
# 运行方式：
# Rscript 02_apply_annotation_plot.R 项目目录 annotation_base.csv 图标题
#
# annotation_base.csv 每行对应一个 cluster，固定包含：
# cluster,base_annotation,prefix_marker
# prefix_marker 从该 cluster 的 Top10 marker 中人工选择。
# 非拟南芥植物只有在一对一或其他高可信唯一映射时使用拟南芥标准 symbol；
# 同源关系不唯一时使用原物种 marker 名称，并在注释依据表中保留全部拟南芥候选。
# 对不重复的基础注释可以留空；脚本只给重复基础注释添加 marker 前缀。

library(Seurat)

# 读取命令行参数，并定位当前脚本所在目录，以便加载同目录的两个辅助脚本。
args <- commandArgs(trailingOnly = TRUE)
project_dir <- args[1]
annotation_file <- args[2]
figure_title <- args[3]

script_argument <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- sub("^--file=", "", script_argument)
script_dir <- dirname(normalizePath(script_path))

source(file.path(script_dir, "resolve_annotation_names.R"))
source(file.path(script_dir, "plot_numbered_umap.R"))

# 读取第一阶段保存的对象和人工填写的基础注释表。
output_dir <- file.path(project_dir, "harmony_analysis")
seurat_object <- readRDS(
  file.path(output_dir, "seurat_harmony_clustered_before_annotation.rds")
)
annotation_table <- read.csv(annotation_file, stringsAsFactors = FALSE)

# 根据基础注释是否重复，生成最终名称。
# 重复名称写成 Marker+ cell type；唯一名称不添加 marker。
annotation_final <- resolve_annotation_names(annotation_table)
write.csv(
  annotation_final,
  file.path(output_dir, "cluster_annotation_final.csv"),
  row.names = FALSE
)

# 将 cluster 到最终注释的映射写入每个细胞的 metadata。
annotation_vector <- setNames(
  annotation_final$final_annotation,
  annotation_final$cluster
)
seurat_object$final_annotation <- annotation_vector[
  as.character(seurat_object$seurat_clusters)
]

# 提取 UMAP 坐标和绘图所需的 metadata，保持一行对应一个细胞。
umap_data <- as.data.frame(Embeddings(seurat_object, reduction = "umap"))
umap_data$cluster <- as.character(seurat_object$seurat_clusters)
umap_data$final_annotation <- seurat_object$final_annotation
umap_data$cell <- rownames(umap_data)

# 保存坐标，便于后续在其他软件中复用或检查。
write.csv(
  umap_data,
  file.path(output_dir, "umap_coordinates_with_annotation.csv"),
  row.names = FALSE
)

# 生成用户指定版式的 PNG、PDF 和 cluster 编号映射表。
plot_numbered_umap(
  umap_data = umap_data,
  title = figure_title,
  png_file = file.path(output_dir, "umap_numbered_annotation.png"),
  pdf_file = file.path(output_dir, "umap_numbered_annotation.pdf"),
  mapping_file = file.path(output_dir, "umap_cluster_number_mapping.csv")
)

# 保存最终注释后的 Seurat 对象。
Idents(seurat_object) <- "final_annotation"
saveRDS(
  seurat_object,
  file.path(output_dir, "seurat_harmony_annotated.rds")
)
