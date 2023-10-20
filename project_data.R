project_data <- function(obj=obj,
                         data_query=df,
                         ref_clusters="sketch_snn_res.2",
                         FSC.A="FSC-A",
                         FSC.H="FSC-H",
                         pred_name="clusters_predicted",
                         chunk_size=1000000){
  
  pkg <- c("Seurat", "tidyr","tidyverse", "dplyr", "BPCells", "readr", "MASS", "pbapply", "data.table")
  invisible(lapply(pkg, library, character.only = TRUE))
  
  data_ref <- as.data.frame(as.data.frame(t(obj@assays$sketch$counts)))
  data_ref$clst <- as.numeric(as.character(obj@meta.data %>% dplyr::filter(!seurat_clusters=="NA") %>% pull(.data[[ref_clusters]])))
  
  lda_model <- MASS::lda(clst ~ ., data=data_ref)
  model_predict <- prediced <- predict(lda_model, data_ref)
  
  if(is.data.frame(data_query) | is.data.table(data_query)){ 
    message("Calculate Ratio")
    data_query=as.data.frame(data_query)
    data_query$ratio <- data_query[,FSC.A]/data_query[,FSC.H]
    data_query$ratio <- scales::rescale(as.numeric(data_query$ratio), to = c(0, 1023))
    data_query <- data_query[is.na(obj$seurat_clusters),]
    
    chunk_size <- chunk_size
    
    if(dim(data_query)[1] < chunk_size){
      chunk_size <- dim(data_query)[1]
    }
    
    n_rows <- dim(data_query)[1]
    message("Calculate Prediction")
    pred <- pblapply(seq(1, n_rows, chunk_size), function(i){
      start <- i
      end_row <- min(i + chunk_size - 1, n_rows)
      data_to_predict <- as.data.frame(data_query[start:end_row,colnames(data_ref)[!colnames(data_ref)=="clst"]])
      prediced <- predict(lda_model, data_to_predict)
      prediced_vector <- as.numeric(as.character(prediced$class))
      return(prediced_vector)
    })
    pred_vector <- unlist(pred)
    #data_query$clusters_predicted <- pred_vector
    obj[[pred_name]] <- "to_predict"
    obj[[pred_name]][is.na(obj$seurat_clusters),] <- as.numeric(as.character(pred_vector))
    obj[[pred_name]][!is.na(obj$seurat_clusters),] <- as.numeric(as.character(data_ref$clst))
    obj[[pred_name]] <- factor(obj@meta.data[,pred_name], levels = sort(unique(as.numeric(obj@meta.data[,pred_name]))))
    return(obj)
    
  }else if(isS4(data_query)){ 
    #data_ref <- obj@meta.data %>% dplyr::filter(!seurat_clusters=="NA") 
    obj_ref <- subset(obj, cells = rownames(data_ref))
    obj_ref <- as.data.frame(t(obj_ref@assays$sketch$counts))
    obj_ref$clst <- as.numeric(as.character(obj@meta.data %>% dplyr::filter(!seurat_clusters=="NA") %>% pull(.data[[ref_clusters]])))
    
    lda_model <- MASS::lda(clst ~ ., data=obj_ref)
    data_to_predict <- as.data.frame(t(as.matrix(obj@assays$FACS$counts)))
    data_to_perdict <- data_to_predict[!rownames(data_to_predict) %in% rownames(obj_ref),]
    
    predicted <- predict(lda_model, data_to_perdict)
    
    obj[[pred_name]] <- "to_predict"
    obj[[pred_name]][!rownames(data_to_predict) %in% rownames(obj_ref),] <- as.numeric(as.character(predicted$class))
    obj[[pred_name]][rownames(data_to_predict) %in% rownames(obj_ref),] <- as.numeric(as.character(data_ref$clst))
    obj[[pred_name]] <-factor(obj@meta.data[,pred_name], levels = sort(unique(as.numeric(obj@meta.data[,pred_name]))))
    return(obj)
  }
}