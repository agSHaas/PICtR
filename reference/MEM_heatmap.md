# Marker Enrichment Modeling (MEM) Heatmap.

Marker Enrichment Modeling (MEM) Heatmap.

## Usage

``` r
MEM_heatmap(
  obj,
  markers = c(),
  cluster_col = "seurat_clusters",
  cols = pals::coolwarm(100),
  heatmap_name = "MEM enrichment score",
  heatmap_column_title = "marker",
  heatmap_row_title = "cluster",
  scale_width = 2.2,
  scale_height = 5
)
```

## Arguments

- obj:

  The Seurat object.

- markers:

  Meta.data columns with features that should be plotted in the heatmap
  and the clustering resolution.

- cluster_col:

  Character string specifying the column that contains the clustering
  solution.

- cols:

  Color palette.

- heatmap_name:

  Title of the heatmap.

- heatmap_column_title:

  Title for the columns of the heatmap.

- heatmap_row_title:

  Title for the rows of the heatmap.

- scale_width:

  Scaling factor for the width of the heatmap in relation to the number
  of columns.

- scale_height:

  Scaling factor for the height of the heatmap in relation to the number
  of rows.

## Value

MEM heat map.
