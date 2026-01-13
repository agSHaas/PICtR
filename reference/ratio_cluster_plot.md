# Ratio cluster plot.

Stacked bar plot of each cluster with the proportion of cells
below/above the determined FSC.A/FSC.H ratio threshold.

## Usage

``` r
ratio_cluster_plot(
  obj,
  clusters = "seurat_clusters",
  ratio = "ratio_anno",
  assay = "FACS"
)
```

## Arguments

- obj:

  The Seurat object.

- clusters:

  The string of the meta.data column with the clustering resolution to
  plot.

- ratio:

  The meta.data column with the classification of cells
  (ratio_high/ratio_low) determined using the FSC.A/FSC.H ratio and the
  determined threshold.

- assay:

  The Seurat assay to use (default FACS).

## Value

None
