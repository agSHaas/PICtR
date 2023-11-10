select_dbt <- function(obj,
                       clusters="seurat_clusters", 
                       ratio="ratio_anno",
                       ratio_high="Ratio_high",
                       assay="FACS",
                       quantile=0.8,
                       selected_clusters="doublet_clusters"){
  
  # pick clusters with a doublet content in the 80th percentile
  dist <- obj@meta.data %>%  group_by(.data[[clusters]], .data[[ratio]]) %>% count(.data[[ratio]]) %>% ungroup() %>% 
    group_by(.data[[clusters]]) %>% mutate(ratio_ratio = n/sum(n)) %>% dplyr::filter(.data[[ratio]]==ratio_high) 
  
  cutoff <- quantile(dist$ratio_ratio, probs = quantile)
  cluster_to_use <- as.vector(pull(dist[dist$ratio_ratio > cutoff,], .data[[clusters]]))
  
  # assign to misc slot in Seurat object as a named vector
  names <- names(obj@misc)
  obj@misc <- append(obj@misc, list(cluster_to_use))
  names(obj@misc) <- append(names, selected_clusters)
  
  # return obj
  return(obj)
}
