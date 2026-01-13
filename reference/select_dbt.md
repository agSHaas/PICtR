# Selection of clusters containing physically interacting cells.

Selects the clusters that belong to a given percentile of the clusters
with the largest physically interacting cells proportion based on FSC
ratio thresholding. By default, the top 20 percent of clusters are
chosen. Consider using
[`ratio_cluster_plot`](https://agshaas.github.io/PICtR/reference/ratio_cluster_plot.md)
to evaluate the choice of percentile cutoff.

## Usage

``` r
select_dbt(
  obj,
  clusters = "seurat_clusters",
  ratio = "ratio_anno",
  ratio_high = "Ratio_high",
  assay = "FACS",
  quantile = 0.8,
  selected_clusters = "doublet_clusters"
)
```

## Arguments

- obj:

  The Seurat object.

- clusters:

  The meta.data column containing the cluster labels (Default:
  seurat_clusters).

- ratio:

  The meta.data column with the classification of cells into ratio_high
  and ratio_low using the FSC ratio (FSC.A/FSC.H) and a thresholding
  method.

- ratio_high:

  The character string that indicates cells with a FSC ratio
  (FSC.A/FSC.H) above the threshold determined with
  [`calculateThreshold`](https://agshaas.github.io/PICtR/reference/calculateThreshold.md)
  within the meta.data column.

- assay:

  The Seurat assay containing FACS data.

- quantile:

  The desired percentile cutoff above which clusters are classified as
  physically interacting cell clusters (Default 0.8, meaning top 20
  percent of clusters are chosen).

- selected_clusters:

  Character vector for the misc slot in the Seurat object that will
  contain the cluster numbers of the selected physically interacting
  cell clusters (Default: "doublet_clusters").

## Value

Seurat object
