# Plot wrapper function.

A wrapper function for plots that are often used during exploratory
analysis within the Seurat framework (<https://satijalab.org/seurat/>).
Includes
[`FeaturePlot`](https://satijalab.org/seurat/reference/FeaturePlot.html),
[`DimPlot`](https://satijalab.org/seurat/reference/DimPlot.html) for
different parameters, and
[`ratio_cluster_plot`](https://agshaas.github.io/PICtR/reference/ratio_cluster_plot.md).

## Usage

``` r
wrapper_for_plots(
  obj = obj,
  feature_plot = TRUE,
  cluster_plot = TRUE,
  meta_list = list("ratio_anno"),
  cluster_handle = "sketch_snn_res",
  feature_plot_colors = pals::parula(1000),
  ratio_plot_color = c(Ratio_low = "dodgerblue2", Ratio_high = "gold2"),
  reduction = "umap",
  alpha = 1,
  raster = F,
  label_size = 3,
  label_box = FALSE,
  assay = "sketch"
)
```

## Arguments

- obj:

  The Seurat object.

- feature_plot:

  Boolean. TRUE indicates that UMAP is colored for all features (see
  [`FeaturePlot`](https://satijalab.org/seurat/reference/FeaturePlot.html))

- cluster_plot:

  Boolean. TRUE indicates that UMAP is colored by the different cluster
  resolutions stored with the cluster_handle.

- meta_list:

  List of meta.data columns to color
  [`DimPlot`](https://satijalab.org/seurat/reference/DimPlot.html) by.
  Can be both numeric to generate
  [`FeaturePlot`](https://satijalab.org/seurat/reference/FeaturePlot.html)
  or class character or factor for
  [`DimPlot`](https://satijalab.org/seurat/reference/DimPlot.html).

- cluster_handle:

  Prefix for the clustering solutions in the meta.data slot.

- feature_plot_colors:

  Color palette for
  [`FeaturePlot`](https://satijalab.org/seurat/reference/FeaturePlot.html).

- ratio_plot_color:

  Color palette for
  [`ratio_cluster_plot`](https://agshaas.github.io/PICtR/reference/ratio_cluster_plot.md).

- reduction:

  Reduction to use for plotting, for example UMAP.

- alpha:

  Alpha value for plotting.

- raster:

  Convert points to raster format. If TRUE plot is rasterized to
  raster.dpi=c(512, 512). Requires ggrastr.

- label_size:

  Size of the labels plotted within the embedding.

- label_box:

  Plot boxes around labels in the color of the cluster.

- assay:

  The Seurat assay (default FACS).

## Value

A list with the requested plots.
