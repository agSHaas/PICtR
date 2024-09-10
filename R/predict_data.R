#' Cluster label prediction
#'
#' Predicts the cluster labels for a reference data set to a query data set using \code{\link[MASS]{lda}}.
#'
#' @param obj The Seurat object.
#' @param data_query A data frame or a Seurat object with the cells whose labels should be predicted.
#' @param ref_clusters The meta.data column with the reference cluster labels.
#' @param FSC.A The name of the column containing the FSC.A scatter parameter.
#' @param FSC.H The name of the column containing the FSC.H scatter parameter.
#' @param pred_name The name of the meta.data column for predicted cluster labels (character vector).
#' @param assay_ref The name of the Seurat assay which was used to calculate the reference cluster labels.
#' @param assay_query The name of the Seurat assay containing cells whose labels should be predicted. Only if the query is provided as a Seurat object.
#' @param chunk_size Chunk size for the prediction progress for verbose output to standard out.
#' @param return_obj Boolean. Add the predicted cluster labels to the Seurat object? Only if the query is provided as a data frame.
#'
#' @return Seurat object or data frame containing the predicted cluster labels.
#'
#' @importFrom MASS lda
#' @importFrom data.table is.data.table
#'
#' @export
predict_data <- project_data <- function(obj=obj,
                         data_query=query,
                         ref_clusters=NULL,
                         FSC.A="FSC.A",
                         FSC.H="FSC.H",
                         pred_name="clusters_predicted",
                         assay_ref=NULL,
                         assay_query=NULL,
                         chunk_size=1000000,
                         return_obj=TRUE){

  # pull clustering with the chosen resolution from sketched cells as training data
  if (is.null(assay_ref)) {
    stop("Please provide assay to pull data from (assay_ref)")
  }
  if (is.null(ref_clusters)) {
    stop("Please provide a column name containing reference cluster labels used as query for projection")
  }
  data_ref <- as.data.frame(t.data.frame(obj[[assay_ref]]$counts))
  data_ref$clst <- as.numeric(as.character(obj@meta.data %>% dplyr::filter(!seurat_clusters=="NA") %>% pull(.data[[ref_clusters]])))

  # train LDA model
  lda_model <- MASS::lda(clst ~ ., data=data_ref)

  # query data as data frame ------------------------------------------------------
  if(is.data.frame(data_query) | data.table::is.data.table(data_query)){
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
    pred <- pbapply::pblapply(seq(1, n_rows, chunk_size), function(i){
      start <- i
      end_row <- min(i + chunk_size - 1, n_rows)
      data_to_predict <- as.data.frame(data_query[start:end_row,colnames(data_ref)[!colnames(data_ref)=="clst"]])
      prediced <- predict(lda_model, data_to_predict)
      prediced_vector <- as.numeric(as.character(prediced$class))
      return(prediced_vector)
    })

    # assign predictions to obj meta.data, keep original clustering for training cells
    pred_vector <- unlist(pred)

    # either return the infomation back to the seurat object or add to the query data
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

  # query data as Seurat object  --------------------------------------------------------
  }else if(isS4(data_query)){

    # pull query data
    if (is.null(assay_query)) {
      stop("Please provide assay to project data for (assay_query)")
    }
    data_to_predict <- as.data.frame(t(as.matrix(obj[[assay_query]]$counts)))
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
