project_data <- function(obj=obj,
                         data_query=df,
                         ref_clusters="",
                         FSC.A="FSC.A",
                         FSC.H="FSC.H",
                         pred_name="clusters_predicted",
                         assay="sketch",
                         chunk_size=1000000, 
                         return_obj=TRUE){
  
  # required packages
  pkg <- c("Seurat", "tidyr","tidyverse", "dplyr", "BPCells", "readr", "MASS", "pbapply", "data.table")
  invisible(lapply(pkg, library, character.only = TRUE))
  
  # pull clustering with the chosen resolution from sketched cells as training data 
  data_ref <- as.data.frame(t(obj[[assay]]$counts))
  data_ref$clst <- as.numeric(as.character(obj@meta.data %>% dplyr::filter(!seurat_clusters=="NA") %>% pull(.data[[ref_clusters]])))
  
  # train LDA model
  lda_model <- MASS::lda(clst ~ ., data=data_ref)
  # predict_ref <- predict(lda_model, data_ref)
  
  # query data as data frame 
  if(is.data.frame(data_query) | is.data.table(data_query)){ 
    data_query=as.data.frame(data_query)
    
    message("Calculate Ratio")
    data_query$ratio <- data_query[,FSC.A]/data_query[,FSC.H]
    data_query$ratio <- scales::rescale(as.numeric(data_query$ratio), to = c(0, 1023))
    
    # remove sketched cells (= training cells)
    if(anyNA(obj$seurat_clusters)){
    data_query <- data_query[is.na(obj$seurat_clusters),]
    }
    
    # chunks 
    chunk_size <- chunk_size
    if(dim(data_query)[1] < chunk_size){
      chunk_size <- dim(data_query)[1]
    }
    
    # predict 
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
    
    # assign predictions to obj meta.data, keep original clustering for training cells
    pred_vector <- unlist(pred)
    #data_query$clusters_predicted <- pred_vector
    
    #Either return the infomation back to the seurat object or add to the query data 
    if(return_obj==TRUE){
      obj[[pred_name]] <- "to_predict"
      obj[[pred_name]][is.na(obj$seurat_clusters),] <- as.numeric(as.character(pred_vector))
      obj[[pred_name]][!is.na(obj$seurat_clusters),] <- as.numeric(as.character(data_ref$clst))
      obj[[pred_name]] <- factor(obj@meta.data[,pred_name], levels = sort(unique(as.numeric(obj@meta.data[,pred_name]))))
      return(obj)
    }else{
      data_query$prediction <- pred_vector
      return(data_query)
    }
    
  # query data as Seurat object  
  }else if(isS4(data_query)){ 
    
    # pull query data 
    data_to_predict <- as.data.frame(t(as.matrix(obj@assays$FACS$counts)))
    data_to_predict <- data_to_predict[!rownames(data_to_predict) %in% rownames(data_ref),]
    
    # prediction
    predicted <- predict(lda_model, data_to_predict)
    
    # add predicted clusters to meta.data, keep original clustering for training cells 
    obj[[pred_name]] <- "to_predict"
    obj[[pred_name]][!colnames(obj) %in% rownames(data_ref),] <- as.numeric(as.character(predicted$class))
    obj[[pred_name]][colnames(obj) %in% rownames(data_ref),] <- as.numeric(as.character(data_ref$clst))
    obj[[pred_name]] <-factor(obj@meta.data[,pred_name], levels = sort(unique(as.numeric(obj@meta.data[,pred_name]))))
    
    return(obj)
  }
}
