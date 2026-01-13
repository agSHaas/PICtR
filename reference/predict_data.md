# Cluster label prediction

Predicts the cluster labels for a reference data set to a query data set
using [`lda`](https://rdrr.io/pkg/MASS/man/lda.html).

## Usage

``` r
predict_data(
  obj = obj,
  data_query = query,
  ref_clusters = NULL,
  FSC.A = "FSC.A",
  FSC.H = "FSC.H",
  pred_name = "clusters_predicted",
  assay_ref = NULL,
  assay_query = NULL,
  chunk_size = 1e+06,
  return_obj = TRUE
)
```

## Arguments

- obj:

  The Seurat object.

- data_query:

  A data frame or a Seurat object with the cells whose labels should be
  predicted.

- ref_clusters:

  The meta.data column with the reference cluster labels.

- FSC.A:

  The name of the column containing the FSC.A scatter parameter.

- FSC.H:

  The name of the column containing the FSC.H scatter parameter.

- pred_name:

  The name of the meta.data column for predicted cluster labels
  (character vector).

- assay_ref:

  The name of the Seurat assay which was used to calculate the reference
  cluster labels.

- assay_query:

  The name of the Seurat assay containing cells whose labels should be
  predicted. Only if the query is provided as a Seurat object.

- chunk_size:

  Chunk size for the prediction progress for verbose output to standard
  out.

- return_obj:

  Boolean. Add the predicted cluster labels to the Seurat object? Only
  if the query is provided as a data frame.

## Value

Seurat object or data frame containing the predicted cluster labels.

## References

Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S.
Fourth Edition. Springer, New York. ISBN 0-387-95457-0

Ripley, B. D. (1996) Pattern Recognition and Neural Networks. Cambridge
University Press.
