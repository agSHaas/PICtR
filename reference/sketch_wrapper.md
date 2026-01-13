# PICtR workflow wrapper

Characterizes cells into ratio_high and ratio_low cells for physically
interacting cell analysis, samples a representative subset of cells
using
[`SketchData`](https://satijalab.org/seurat/reference/SketchData.html)
and runs the standard Seurat analysis workflow.

## Usage

``` r
sketch_wrapper(
  channel = channel,
  meta_data = NULL,
  assay = "FACS",
  FSC.A = "FSC.A",
  FSC.H = "FSC.H",
  n_sketch_cells = 50000,
  n_components = "all",
  clst_algorithm = 1,
  resolution = c(0.5, 1, 2, 3, 4),
  min_clst_size = 100,
  meta.k = "auto",
  obj_name = "obj_sketched_non_projected",
  ratio = TRUE,
  thresholding_method = "otsu",
  hist_breaks = 2000,
  group_by = NULL,
  verbose = TRUE,
  BPcell_dir = NULL,
  overwrite = F,
  working_dir = getwd(),
  seed = 42
)
```

## Arguments

- channel:

  A data frame with the dimensionality cells x flow cytometry
  parameters.

- meta_data:

  A data frame with meta_data for every cell (Default=NULL).

- assay:

  A character string with the name of the assay to be created
  (Default="FACS").

- FSC.A:

  The name (string) of the column containing the FSC.A scatter parameter
  (Default="FACS.A").

- FSC.H:

  The name (string) of the column containing the FSC.H scatter parameter
  (Default="FACS.H").

- n_sketch_cells:

  The number of cells to be subsampled by
  [`SketchData`](https://satijalab.org/seurat/reference/SketchData.html)
  (Default=50000).

- n_components:

  The number of components computed and used during analysis. "all" or
  given as an integer (Default="all").

- clst_algorithm:

  Algorithm used for clustering. For Seurat's
  [`FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html):
  1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel
  refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires
  the leidenalg python. Alternatively, "hdbscan" for Hierarchical
  Density-Based Spatial Clustering of Applications with Noise (HDBSCAN)
  on the UMAP embedding, see
  [`hdbscan`](https://rdrr.io/pkg/dbscan/man/hdbscan.html)), "flowMeans"
  for non-parametric clustering and segmented-regression-based change
  point detection, see
  [`flowMeans`](https://rdrr.io/pkg/flowMeans/man/flowMeans.html), or
  "flowSOM" for self-organizing maps, see
  [`run.flowsom`](http://immunedynamics.io/spectre/reference/run.flowsom.md),
  (Default: Louvain, 1)

- resolution:

  A numerical vector with the desired resolutions for clustering
  (Default: c(0.5,1,2,3,4)). Only used for Louvain, SLM, and Leiden
  clustering.

- min_clst_size:

  minimum cluster size for HDBSCAN. See details at
  [`hdbscan`](https://rdrr.io/pkg/dbscan/man/hdbscan.html) (Default:
  100).

- meta.k:

  meta k for FlowSOM clustering, see
  [`run.flowsom`](http://immunedynamics.io/spectre/reference/run.flowsom.md)
  (Default: "auto").

- obj_name:

  The name used for storing the Seurat object (character).

- ratio:

  Boolean variable. If TRUE the FSC ratio (FSC.A/FSC.H) will be
  calculated.

- thresholding_method:

  Method for thresholding. One of "otsu", "triangle", "kmeans", or any
  method from the autothresholdr package. For details see
  ?calculateThreshold() (Default = "otsu").

- hist_breaks:

  Number of histogram breaks for methods that rely on histograms. Should
  be \>100 (Default: 2000).

- group_by:

  Optional grouping parameter to calculate the FSC ratio (FSC.A/FSC.H)
  thresholding method for.

- verbose:

  Verbosity (Boolean).

- BPcell_dir:

  Optional directory with the counts matrix for
  [`open_matrix_dir`](https://bnprks.github.io/BPCells/reference/matrix_io.html).

- overwrite:

  Overwrite existing BPCells directory? (Default: FALSE)

- working_dir:

  Directory path to be used as working directory (character string).

- seed:

  Seed for reproducibility (Default: 42).

## Value

Seurat object.

## References

Hao et al. Dictionary learning for integrative, multimodal and scalable
single-cell analysis. Nature Biotechnology (2023). doi:
<https://doi.org/10.1038/s41587-023-01767-y>.

Van Gassen S et al. (2015) FlowSOM: Using self-organizing maps for
visualization and interpretation of cytometry data. Cytom Part J Int Soc
Anal Cytol 87: 636-645.
<https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625>.

Aghaeepour, N., Nikolic, R., Hoos, H.H. and Brinkman, R.R. (2011), Rapid
cell population identification in flow cytometry data. Cytometry, 79A:
6-13. <https://doi.org/10.1002/cyto.a.21007>.

Campello, R.J.G.B., Moulavi, D., Sander, J. (2013). Density-Based
Clustering Based on Hierarchical Density Estimates. In: Pei, J., Tseng,
V.S., Cao, L., Motoda, H., Xu, G. (eds) Advances in Knowledge Discovery
and Data Mining. PAKDD 2013. Lecture Notes in Computer Science(), vol
7819. Springer, Berlin, Heidelberg.
<https://doi.org/10.1007/978-3-642-37456-2_14>

## Examples

``` r
obj <- sketch_wrapper(channel = demo_lcmv,
                      meta_data = demo_lcmv,
                      n_sketch_cells = 5000,
                      clst_algorithm = 1,
                      ratio = TRUE,
                      thresholding_method = "otsu", 
                      working_dir = tempdir())
#> Calculation of Ratio
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Matrix compression performs poorly with non-integers.
#> • Consider calling convert_matrix_type if a compressed integer matrix is intended.
#> This message is displayed once every 8 hours.
#> Your newly generated object will be saved under: /tmp/RtmpbMnwdx/obj_sketched_non_projected.rds
#> Sketching started...
#> Finding variable features for layer counts
#> Warning: Chernobyl! trL>n 15
#> Warning: Chernobyl! trL>n 15
#> Warning: NaNs produced
#> Calcuating Leverage Score
#> Attempting to cast layer counts to dgCMatrix
#> Attempting to cast layer data to dgCMatrix
#> Sketching is done
#> Finding variable features for layer counts
#> Warning: Chernobyl! trL>n 15
#> Warning: Chernobyl! trL>n 15
#> Warning: NaNs produced
#> Centering and scaling data matrix
#> Warning: Requested number is larger than the number of available items (15). Setting to 15.
#> Warning: Requested number is larger than the number of available items (15). Setting to 15.
#> Warning: Requested number is larger than the number of available items (15). Setting to 15.
#> Warning: Requested number is larger than the number of available items (15). Setting to 15.
#> Warning: Requested number is larger than the number of available items (15). Setting to 15.
#> PC_ 1 
#> Positive:  CD4, Ly6G, CD11b, CD90-1, CD8, CD11c, SSC.H, MHCII 
#> Negative:  FSC.A, ratio, SSC.A, FSC.H, CD45-2, CD3, CD19 
#> PC_ 2 
#> Positive:  CD4, CD45-2, CD3, CD90-1, CD8, MHCII, CD19, ratio 
#> Negative:  CD11b, SSC.H, Ly6G, SSC.A, FSC.H, CD11c, FSC.A 
#> PC_ 3 
#> Positive:  CD4, CD90-1, CD3, SSC.H, CD11b, SSC.A, FSC.H, CD8 
#> Negative:  MHCII, CD19, ratio, CD45-2, FSC.A, CD11c, Ly6G 
#> PC_ 4 
#> Positive:  CD19, CD90-1, CD4, Ly6G, ratio, FSC.A, SSC.A, MHCII 
#> Negative:  CD11c, CD8, CD45-2, CD11b, CD3, SSC.H, FSC.H 
#> PC_ 5 
#> Positive:  CD8, Ly6G, CD19, CD3, FSC.H, FSC.A, ratio, SSC.A 
#> Negative:  CD11c, MHCII, CD4, CD90-1, CD45-2, SSC.H, CD11b 
#> Computing nearest neighbor graph
#> Computing SNN
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 5000
#> Number of edges: 164643
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.9193
#> Number of communities: 14
#> Elapsed time: 0 seconds
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 5000
#> Number of edges: 164643
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.8813
#> Number of communities: 19
#> Elapsed time: 0 seconds
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 5000
#> Number of edges: 164643
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.8238
#> Number of communities: 24
#> Elapsed time: 0 seconds
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 5000
#> Number of edges: 164643
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.7740
#> Number of communities: 32
#> Elapsed time: 0 seconds
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 5000
#> Number of edges: 164643
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.7367
#> Number of communities: 38
#> Elapsed time: 0 seconds
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
#> UMAP will return its model
#> 17:30:41 UMAP embedding parameters a = 0.9922 b = 1.112
#> 17:30:41 Read 5000 rows and found 15 numeric columns
#> 17:30:41 Using Annoy for neighbor search, n_neighbors = 30
#> 17:30:41 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> *
#> |
#> 17:30:41 Writing NN index file to temp file /tmp/RtmpbMnwdx/file9b1a4213f4d7
#> 17:30:41 Searching Annoy index using 1 thread, search_k = 3000
#> 17:30:42 Annoy recall = 100%
#> 17:30:43 Commencing smooth kNN distance calibration using 1 thread
#>  with target n_neighbors = 30
#> 17:30:43 Initializing from normalized Laplacian + noise (using RSpectra)
#> 17:30:43 Commencing optimization for 500 epochs, with 202814 positive edges
#> 17:30:43 Using rng type: pcg
#> 17:30:48 Optimization finished
#> The object will be updated and saved
```
