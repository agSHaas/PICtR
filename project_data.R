project_data <- function(obj,
                         ref_clusters="sketch_snn_res.2",
                         return=c("obj", "df")){
  pkg <- c("Seurat", "tidyr","tidyverse", "dplyr", "BPCells", "readr", "MASS")
  invisible(lapply(pkg, library, character.only = TRUE))
  data_ref <- obj@meta.data %>% filter(!seurat_clusters=="NA")
  obj_ref <- subset(obj, cells = rownames(data_ref))
  obj_ref <- as.data.frame(t(obj_ref@assays$sketch$counts))
  obj_ref$clst <- as.vector(data_ref[,ref_clusters])
  
  
  obj_query <- subset(obj, cells = rownames(obj@meta.data %>% filter(is.na(seurat_clusters))))
  obj_query <- as.data.frame(t(as.matrix(obj_query@assays$FACS$counts)))
  
  lda_model <- MASS::lda(clst ~ ., data=obj_ref)
  
  data_to_perdict <- obj_query[,colnames(obj_ref)[-"clst"]]
  
  prediced <- predict(lda_model, data_to_perdict)
  
  if(return=="obj"){
    obj$cluster_predicted <- factor(prediced$class, levels = sort(as.numeric(as.vector(unique(prediced$class)))))
    return(obj)
  }else if(retrun=="df"){
    obj_query$cluster_predicted <- factor(prediced$class, levels = sort(as.numeric(as.vector(unique(prediced$class)))))
    obj_ref$status <- "ref"
    obj_query[,colnames(obj_ref)[-"clst"]]$status <- "query"
    obj_df <- rbind(obj_ref, obj_query[,colnames(obj_ref)[-"clst"]]$status)
    return(obj_df)
  }
  

}
