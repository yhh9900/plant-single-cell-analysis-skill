# 植物 marker 与拟南芥基因名匹配

## 基本原则

拟南芥基因名用于跨物种解释和细胞类型注释，不用于替换原始表达矩阵的 feature 名称。任何输出都同时保留原物种 ID、原物种 symbol、拟南芥 locus ID 和拟南芥标准 symbol。

## 拟南芥数据

从本次计数所用的同一版本 GTF/GFF 或配套 gene-info 表建立映射。推荐至少保存：

- `source_gene_id`：矩阵中的原始 feature 名称；
- `source_gene_name`：参考注释中的基因名；
- `arabidopsis_gene_id`：TAIR locus ID，例如 `AT1G01010`；
- `arabidopsis_gene_name`：标准 symbol；
- `orthology_type`：拟南芥自身数据可写 `self`；
- `mapping_source`：注释版本或文件来源。

不要只保存 symbol，因为同义词、旧名称和缺失 symbol 会造成 marker 无法追溯。

## 其他植物

使用与当前参考蛋白/基因版本一致的同源表。若项目已有 OrthoFinder 结果，优先从对应 orthogroup 和一对一 ortholog 表生成映射；不要把不同版本的基因 ID 混在一起。

同源关系按以下方式处理：

- 一对一：可直接使用拟南芥 symbol 辅助注释。
- 一对多：保留所有拟南芥候选，不自动挑选最熟悉的名称。
- 多对一：多个原物种基因分别保留；不能先求和再计算 marker。
- 多对多：完整保存候选关系，并结合表达特异性、蛋白同源证据和已发表植物 marker 判断。
- 多倍体：保留各亚基因组 homeolog，尤其不能把棉花 homeolog 静默折叠成一个拟南芥基因。

## Top10 marker 表

Top10 的排名基于原物种差异表达结果。先选出每个 cluster 的 10 个原物种 marker，再连接拟南芥映射表；一对多映射会使输出行数超过 10，但仍只代表 10 个原物种 marker。

推荐输出列：

`cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2, source_gene_name, arabidopsis_gene_id, arabidopsis_gene_name, orthology_type, mapping_source`

## UMAP 注释前缀

重复细胞类型需要 marker 前缀时：

- 唯一且高可信的拟南芥匹配：可使用拟南芥标准 symbol，便于跨物种统一展示。
- 拟南芥匹配不唯一：使用原物种 marker 名称，不把某一个候选伪装成唯一对应关系。
- 无 symbol：保留稳定 locus ID，不自行编造缩写。
