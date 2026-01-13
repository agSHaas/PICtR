# Split plot wrapper.

Dimensional reduction (UMAP) plot split by a given parameter. Per
default the returned split plots are rasterized using
[`geom_point_rast`](https://rdrr.io/pkg/ggrastr/man/geom_point_rast.html).
Requires ggrastr.

## Usage

``` r
split_plot_sketch(obj, group_by = "seurat_clusters", split_by = "ratio_anno")
```

## Arguments

- obj:

  The Seurat object.

- group_by:

  Parameter to group the plot by.

- split_by:

  Parameter to split the plot by.

## Value

ggplot object
