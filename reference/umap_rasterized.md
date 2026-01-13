# Rasterized UMAP

Wrapper function to plot a dimensional reduction plot colored by
clusters.

## Usage

``` r
umap_rasterized(
  data = obj,
  group.by = "seurat_clusters",
  raster.dpi = 500,
  label = TRUE,
  cols = pals::tol.rainbow(70),
  umap1 = "umap_1",
  umap2 = "umap_2",
  reduction = "umap",
  raster = TRUE
)
```

## Arguments

- data:

  Seurat object.

- group.by:

  meta.data column to group the dimensional reduction plot by, for
  example a clustering solution.

- raster.dpi:

  Pixel resolution (numeric; default 500).

- label:

  Boolean. Plot labels?

- cols:

  Color palette.

- umap1:

  First UMAP dimension as found in the meta.data.

- umap2:

  Second UMAP dimension as found in the meta.data.

- reduction:

  Which reduction to use.

- raster:

  Boolean. Rasterize the plot?

## Value

Plot
