# 功能：把原物种 marker 表连接到拟南芥 locus ID 和标准基因名。
#
# 运行方式：
# Rscript map_markers_to_arabidopsis.R markers_all.csv gene_to_arabidopsis_mapping.csv 输出目录
#
# gene_to_arabidopsis_mapping.csv 固定包含：
# source_gene_id,source_gene_name,arabidopsis_gene_id,arabidopsis_gene_name,orthology_type,mapping_source
#
# source_gene_id 必须与 Seurat marker 表中的 gene 列完全一致。
# 一对多或多对多关系在映射表中使用多行表示，连接后也会保留多行。
# 本脚本不修改 Seurat 对象的行名，也不合并任何基因的表达量。

library(dplyr)

# 读取命令行参数和两个输入表。
args <- commandArgs(trailingOnly = TRUE)
markers_file <- args[1]
mapping_file <- args[2]
output_dir <- args[3]

markers_all <- read.csv(markers_file, stringsAsFactors = FALSE)
gene_mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)

# 保存原始排序位置，连接一对多映射后仍按原 marker 顺序输出。
markers_all$marker_order <- seq_len(nrow(markers_all))

# Top10 排名先在原物种 marker 上完成。
# 一对多映射发生在 Top10 选择之后，因此不会让拥有多个拟南芥候选的基因获得额外排名优势。
markers_top10 <- markers_all %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj) %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

# 左连接保留全部原物种 marker；没有拟南芥对应关系的 marker 仍保留在输出中。
markers_all_mapped <- merge(
  markers_all,
  gene_mapping,
  by.x = "gene",
  by.y = "source_gene_id",
  all.x = TRUE,
  sort = FALSE
)
markers_all_mapped <- markers_all_mapped[
  order(markers_all_mapped$marker_order),
]

markers_top10_mapped <- merge(
  markers_top10,
  gene_mapping,
  by.x = "gene",
  by.y = "source_gene_id",
  all.x = TRUE,
  sort = FALSE
)
markers_top10_mapped <- markers_top10_mapped[
  order(markers_top10_mapped$cluster, markers_top10_mapped$marker_order),
]

# 生成便于阅读的拟南芥 marker 标签。
# 有标准 symbol 时显示“SYMBOL (AT locus)”，只有 locus ID 时保留稳定 ID。
markers_all_mapped$arabidopsis_marker <- ifelse(
  is.na(markers_all_mapped$arabidopsis_gene_name) |
    markers_all_mapped$arabidopsis_gene_name == "",
  markers_all_mapped$arabidopsis_gene_id,
  paste0(
    markers_all_mapped$arabidopsis_gene_name,
    " (",
    markers_all_mapped$arabidopsis_gene_id,
    ")"
  )
)
markers_top10_mapped$arabidopsis_marker <- ifelse(
  is.na(markers_top10_mapped$arabidopsis_gene_name) |
    markers_top10_mapped$arabidopsis_gene_name == "",
  markers_top10_mapped$arabidopsis_gene_id,
  paste0(
    markers_top10_mapped$arabidopsis_gene_name,
    " (",
    markers_top10_mapped$arabidopsis_gene_id,
    ")"
  )
)

# 统计每个 cluster 的 Top10 原物种 marker 中有多少能匹配到拟南芥。
mapping_summary <- markers_top10_mapped %>%
  group_by(cluster) %>%
  summarise(
    source_top10_markers = n_distinct(gene),
    mapped_source_markers = n_distinct(gene[!is.na(arabidopsis_gene_id)]),
    .groups = "drop"
  )

# 输出完整 marker、Top10 marker 和每个 cluster 的映射覆盖率。
write.csv(
  markers_all_mapped,
  file.path(output_dir, "markers_all_with_arabidopsis.csv"),
  row.names = FALSE
)
write.csv(
  markers_top10_mapped,
  file.path(output_dir, "markers_top10_with_arabidopsis.csv"),
  row.names = FALSE
)
write.csv(
  mapping_summary,
  file.path(output_dir, "arabidopsis_mapping_summary_by_cluster.csv"),
  row.names = FALSE
)
