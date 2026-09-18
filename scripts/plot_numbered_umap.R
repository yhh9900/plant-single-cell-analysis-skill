# 功能：绘制“左侧 UMAP 点图 + 右侧编号注释区”的组合图。
# 图形结构按照用户给出的示例设计：
#   - 左侧只显示彩色细胞点和 cluster 编号，不在点云上放长注释。
#   - 右侧使用带编号的彩色圆点，并在圆点右边显示最终注释。
#   - 上方显示居中的项目标题。
#
# 输入 umap_data 每行代表一个细胞，固定包含：
#   UMAP_1、UMAP_2、cluster、final_annotation。
# cluster 可以是 Seurat 的 0、1、2……，图中会另外转换为 1、2、3……展示编号。
# 原始 cluster ID 不会被修改，并会写入编号映射表。

library(ggplot2)
library(patchwork)

plot_numbered_umap <- function(
  umap_data,
  title,
  png_file,
  pdf_file,
  mapping_file,
  width = 10,
  height = 7,
  point_size = 0.25,
  point_alpha = 0.85
) {
  # 提取原始 cluster ID，并按数值顺序排列。
  # 例如 Seurat cluster 0、1、2 会对应图中的展示编号 1、2、3。
  cluster_levels <- unique(as.character(umap_data$cluster))
  cluster_levels <- cluster_levels[order(as.numeric(cluster_levels))]

  # 把 cluster 转为有固定顺序的因子，保证点图、颜色和右侧注释顺序一致。
  umap_data$cluster <- factor(umap_data$cluster, levels = cluster_levels)

  # 为每个原始 cluster 建立连续的 1-based 展示编号。
  cluster_mapping <- data.frame(
    original_cluster = cluster_levels,
    display_number = seq_along(cluster_levels),
    stringsAsFactors = FALSE
  )

  # 每个 cluster 的最终注释在所有细胞中相同，因此取第一条即可。
  annotation_by_cluster <- tapply(
    as.character(umap_data$final_annotation),
    umap_data$cluster,
    function(x) x[1]
  )
  cluster_mapping$final_annotation <-
    annotation_by_cluster[cluster_mapping$original_cluster]

  # 使用色相分布均匀的离散颜色，并把颜色名称绑定到原始 cluster ID。
  cluster_colors <- grDevices::hcl.colors(
    length(cluster_levels),
    palette = "Dynamic"
  )
  names(cluster_colors) <- cluster_levels
  cluster_mapping$color <- cluster_colors[cluster_mapping$original_cluster]

  # 将展示编号映射到每一个细胞，后面用于计算编号位置。
  umap_data$display_number <- cluster_mapping$display_number[
    match(as.character(umap_data$cluster), cluster_mapping$original_cluster)
  ]

  # 使用每个 cluster 的 UMAP 坐标中位数作为编号位置。
  # 中位数比均值更不容易被少量离群细胞拉偏。
  cluster_centers <- aggregate(
    cbind(UMAP_1, UMAP_2) ~ cluster,
    data = umap_data,
    FUN = median
  )
  cluster_centers$display_number <- cluster_mapping$display_number[
    match(as.character(cluster_centers$cluster), cluster_mapping$original_cluster)
  ]

  # 左侧 UMAP：白底、经典坐标轴、小点、cluster 中心黑色编号。
  umap_panel <- ggplot(
    umap_data,
    aes(x = UMAP_1, y = UMAP_2, color = cluster)
  ) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_text(
      data = cluster_centers,
      aes(label = display_number),
      color = "black",
      fontface = "bold",
      size = 3.5,
      show.legend = FALSE
    ) +
    scale_color_manual(values = cluster_colors, guide = "none") +
    labs(x = "UMAP_1", y = "UMAP_2") +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.line = element_line(color = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )

  # 右侧注释区从上到下按展示编号排列。
  # y 值反向设置，使编号 1 出现在最上方。
  legend_data <- cluster_mapping
  legend_data$y <- rev(seq_len(nrow(legend_data)))
  legend_data$x <- 0
  legend_data$original_cluster <- factor(
    legend_data$original_cluster,
    levels = cluster_levels
  )

  # 彩色圆点用 shape 21 绘制，从而同时控制内部填色和黑色边框。
  # 随后在圆点中心叠加编号，在右侧写最终注释。
  annotation_panel <- ggplot(legend_data, aes(x = x, y = y)) +
    geom_point(
      aes(fill = original_cluster),
      shape = 21,
      size = 5.6,
      stroke = 0.6,
      color = "black"
    ) +
    geom_text(aes(label = display_number), size = 2.8, color = "black") +
    geom_text(
      aes(x = 0.18, label = final_annotation),
      hjust = 0,
      size = 3.4,
      color = "black"
    ) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    coord_cartesian(
      xlim = c(-0.12, 1.8),
      ylim = c(0.35, nrow(legend_data) + 0.65),
      clip = "off"
    ) +
    theme_void() +
    theme(plot.margin = margin(8, 12, 8, 8))

  # 使用 patchwork 横向拼接两部分，并在整张图上方添加居中标题。
  combined_plot <- umap_panel + annotation_panel +
    plot_layout(widths = c(1.25, 1)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 16
        )
      )
    )

  # 同时保存高分辨率 PNG、矢量 PDF 和 cluster 编号映射表。
  ggsave(
    png_file,
    combined_plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    pdf_file,
    combined_plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
  write.csv(cluster_mapping, mapping_file, row.names = FALSE)

  # 返回图对象和映射表，方便在 R 会话中继续修改或检查。
  list(plot = combined_plot, mapping = cluster_mapping)
}
