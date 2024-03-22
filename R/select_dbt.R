#' Selection of clusters containing physically interacting cells.
#'
#' Selects the clusters that belong to a given percentile of the clusters with the largest physically interacting cells proportion based on
#'
#' @param obj The Seurat object.
#' @param clusters The meta.data column with the clustering solution.
#' @param ratio The meta.data column with the classification of cells based on the Otsu threshold of the FSC.A/FSC.H ratio.
#' @param ratio_high The character string that indicates cells with a FSC.A/FSC.H ratio above the Otsu threshold within the meta.data column with the classification of cells based on the FSC.A/FSC.H ratio.
#' @param assay The Seurat assay.
#' @param quantile The desired percentile cutoff above which clusters are classified as physically interacting cell clusters.
#' @param selected_clusters Character vector for the misc slot in the Seurat object that will contain the selected physically interacting cell clusters.
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
  dist <- obj@meta.data %>%  group_by(.data[[clusters]], .data[[ratio]]) %>% count(.data[[ratio]]) %>% ungroup() %>%
    group_by(.data[[clusters]]) %>% mutate(ratio_ratio = n/sum(n)) %>% dplyr::filter(.data[[ratio]]==ratio_high)

  cutoff <- quantile(dist$ratio_ratio, probs = quantile)
  cluster_to_use <- as.vector(pull(dist[dist$ratio_ratio > cutoff,], .data[[clusters]]))

  # assign to misc slot in Seurat object as a named vector
  names <- names(obj@misc)
  obj@misc <- append(obj@misc, list(cluster_to_use))
  names(obj@misc) <- append(names, paste0(selected_clusters, "_q", quantile))

  # return obj
  return(obj)
}
