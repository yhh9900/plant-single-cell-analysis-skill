---
name: single-cell-harmony-annotation
description: Process plant single-cell or single-nucleus RNA-seq from FASTQ through count generation, plant-aware QC, Seurat/Harmony batch correction, clustering, Arabidopsis ortholog/name matching, Top10-marker review, marker-based annotation, and publication-ready numbered UMAP plots. Use for end-to-end plant sc/snRNA workflows; do not use for animal or bulk RNA-seq.
metadata:
  short-description: FASTQ to Harmony clustering and annotated UMAP
---

# Single-cell Harmony Annotation

Build a reproducible plant sc/snRNA-seq workflow from FASTQ to an annotated UMAP. Use Seurat for downstream analysis and Harmony for batch correction unless the user explicitly selects another implementation.

## Resolve inputs before expensive processing

Confirm or infer from reliable metadata:

- assay type and chemistry, including which read contains barcode, UMI, and transcript;
- single-cell versus single-nucleus data;
- sample names, biological groups, and the technical batch variable for Harmony;
- reference genome, annotation release, and expected gene identifiers;
- plant species and the table used to map native genes to Arabidopsis thaliana locus IDs and standard gene names;
- whether a validated count matrix already exists. If it does, resume from the matrix instead of recounting FASTQ.

Do not guess chemistry-specific barcode/UMI settings from filenames alone. Preserve technical reads during FASTQ conversion and inspect read lengths when metadata is incomplete.

## Required workflow

Read [references/workflow.md](references/workflow.md) before running the analysis. For Arabidopsis naming or cross-species plant annotation, also read [references/plant-arabidopsis-mapping.md](references/plant-arabidopsis-mapping.md).

1. Validate FASTQ pairs, read structure, checksums when available, and per-sample metadata.
2. Generate a gene-by-cell count matrix with the chemistry-appropriate counter. For single-nucleus data, count intronic reads when supported by the assay and reference.
3. Perform sample-aware plant QC, including separate mitochondrial and chloroplast fractions, then normalize, select variable features, scale, and run PCA. Record all thresholds and retained cell counts.
4. Run Harmony on a technical batch field. Build neighbors, clusters, and UMAP from the Harmony reduction, not from uncorrected PCA.
5. Compute markers for every cluster before assigning cell types. Export all markers and a ranked Top10 marker table per cluster.
6. Keep native marker IDs and add Arabidopsis locus ID and standard gene-name columns using the project ortholog table. Do not replace or collapse the expression-matrix row names.
7. Assign a base biological annotation from the native and Arabidopsis-matched Top10 markers, the full marker table when needed, and plant references. Keep an evidence table linking every cluster to supporting and conflicting markers.
8. Resolve final annotation names with the duplicate-label rule below.
9. Render the numbered UMAP using [references/umap-style.md](references/umap-style.md) and `scripts/plot_numbered_umap.R`.

## Annotation naming rule

Assign base annotations first, then evaluate duplicates across clusters.

- If a base annotation occurs once, keep only the cell-type name. Do not add a marker prefix.
- If a base annotation occurs in multiple clusters, prefix each repeated label with one discriminating marker from that cluster, formatted exactly as `Marker+ cell type`.
- Prefer a marker from that cluster's Top10 list that distinguishes it from the other clusters sharing the same base annotation.
- For non-Arabidopsis species, use the Arabidopsis standard symbol as the displayed prefix only when the mapping is unique and well supported. For ambiguous one-to-many mappings, retain the native marker name in the figure and document all Arabidopsis candidates in the evidence table.
- If the Top10 list has no defensible marker, inspect the full ranked marker list and document the choice. Never invent a marker.

Use `scripts/resolve_annotation_names.R` after expert review of base annotations and marker candidates.

## Code style

Write direct research scripts around a clearly documented input contract. Do not add defensive-programming layers, broad exception handling, repeated existence checks, or automatic fallback branches. Add detailed Chinese comments explaining the purpose, inputs, outputs, assumptions, and adjustable parameters of each analysis block.

## Output contract

Produce, at minimum:

- FASTQ/counting QC report and sample manifest;
- raw and filtered count matrices or links to the counter output;
- a processed Seurat object with `pca`, `harmony`, and `umap` reductions;
- QC thresholds and cell counts by sample before and after filtering;
- `markers_all.csv` and `markers_top10_by_cluster.csv`;
- `gene_to_arabidopsis_mapping.csv`, `markers_all_with_arabidopsis.csv`, and `markers_top10_with_arabidopsis.csv`;
- `cluster_annotation_evidence.csv` and `cluster_annotation_final.csv`;
- UMAP coordinates and the original-cluster-to-display-number mapping;
- numbered annotated UMAP as high-resolution PNG and vector PDF;
- parameters, software versions, reference versions, and `sessionInfo()`.

Preserve original cluster IDs in tables and objects even when the figure uses display numbers `1..N`.
