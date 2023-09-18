select_dbt <- function(obj,
                       clusters="seurat_clusters", 
                       ratio="ratio_anno",
                       assay="FACS"){
  dist <- obj@meta.data %>%  group_by(.data[[clusters]], .data[[ratio]]) %>% count(.data[[ratio]]) %>% ungroup() %>% 
    group_by(.data[[clusters]]) %>% mutate(ratio_ratio = n/sum(n)) %>% dplyr::filter(.data[[ratio]]=="Ratio_high") 
  
  cutoff <- quantile(dist$ratio_ratio, probs = 80/100)
  cluster_to_use <- as.vector(pull(dist[dist$ratio_ratio > cutoff,], .data[[clusters]]))
  return(cluster_to_use)
}