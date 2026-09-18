# 功能：根据基础细胞类型是否重复，生成最终注释名称。
#
# 输入 annotation_table 每行代表一个 cluster，固定包含三列：
#   cluster          原始 cluster ID。
#   base_annotation  根据 Top10 marker 判断出的基础细胞类型。
#   prefix_marker    从该 cluster 的 Top10 marker 中人工确认的区分 marker。
#                    只有基础注释重复的 cluster 会使用这一列。
#
# 命名规则：
#   1. 基础注释只出现一次：最终名称仍是 base_annotation。
#   2. 基础注释出现多次：最终名称写成 Marker+ cell type。
#
# 这个函数不自动判断 marker 的生物学意义。prefix_marker 应在查看 Top10 marker、
# 完整 marker 表和参考资料后填写，避免用纯算法替代人工注释判断。

resolve_annotation_names <- function(annotation_table) {
  # 统计每个基础细胞类型一共被多少个 cluster 使用。
  annotation_frequency <- table(annotation_table$base_annotation)

  # 把统计结果映射回每一行；大于 1 表示该基础注释发生了重复。
  annotation_table$is_duplicated <-
    annotation_frequency[annotation_table$base_annotation] > 1

  # 先把基础注释复制为最终注释，保证不重复的名称不添加 marker。
  annotation_table$final_annotation <- annotation_table$base_annotation

  # 对重复注释添加一个经过人工确认的区分 marker。
  duplicated_rows <- annotation_table$is_duplicated
  annotation_table$final_annotation[duplicated_rows] <- paste0(
    annotation_table$prefix_marker[duplicated_rows],
    "+ ",
    annotation_table$base_annotation[duplicated_rows]
  )

  # 返回完整映射表，后续用于写入 Seurat metadata 和绘制 UMAP。
  annotation_table
}
