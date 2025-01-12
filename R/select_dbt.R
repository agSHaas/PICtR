#' Selection of clusters containing physically interacting cells.
#'
#' Selects the clusters that belong to a given percentile of the clusters with the largest physically interacting cells proportion based on FSC ratio thresholding. 
#' By default, the top 20 percent of clusters are chosen. Consider using \code{\link{ratio_cluster_plot}} to evaluate the choice of percentile cutoff.
#'
#' @param obj The Seurat object.
#' @param clusters The meta.data column containing the cluster labels (Default: seurat_clusters).
#' @param ratio The meta.data column with the classification of cells into ratio_high and ratio_low using the FSC ratio (FSC.A/FSC.H) and a thresholding method.
#' @param ratio_high The character string that indicates cells with a FSC ratio (FSC.A/FSC.H) above the threshold determined with \code{\link{calculateThreshold}} within the meta.data column.
#' @param assay The Seurat assay containing FACS data.
#' @param quantile The desired percentile cutoff above which clusters are classified as physically interacting cell clusters (Default 0.8, meaning top 20 percent of clusters are chosen).
#' @param selected_clusters Character vector for the misc slot in the Seurat object that will contain the cluster numbers of the selected physically interacting cell clusters (Default: "doublet_clusters").
#'
#' @return Seurat object
#'
#' @export
select_dbt <- function(obj,
                       clusters="seurat_clusters",
                       ratio="ratio_anno",
                       ratio_high="Ratio_high",
                       assay="FACS",
                       quantile=0.8,
                       selected_clusters="doublet_clusters"){

  # pick clusters with a doublet content in the given percentile
  message("Choosing top 20 % of clusters...")
  dist <- obj@meta.data %>%  group_by(.data[[clusters]], .data[[ratio]]) %>% count(.data[[ratio]]) %>% ungroup() %>%
    group_by(.data[[clusters]]) %>% mutate(ratio_ratio = n/sum(n)) %>% dplyr::filter(.data[[ratio]]==ratio_high)

  cutoff <- quantile(dist$ratio_ratio, probs = quantile)
  cluster_to_use <- as.vector(pull(dist[dist$ratio_ratio > cutoff,], .data[[clusters]]))

  # assign to misc slot in Seurat object as a named vector
  message(paste0("Chosen clusters: ", paste0(cluster_to_use, collapse = ",")))
  names <- names(obj@misc)
  obj@misc <- append(obj@misc, list(cluster_to_use))
  names(obj@misc) <- append(names, paste0(selected_clusters, "_q", quantile))

  # return obj
  return(obj)
}
