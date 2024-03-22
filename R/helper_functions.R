#' Save Seurat object
#'
#' Wrapper function for saving Seurat objects
#'
#' @param object The Seurat Object.
#' @param path Directory path.
#' @param file File name.
#'
#' @return None
#'
#' @export
save.obj <- function(object, path=output_obj_path, file){
  if(is.null(file)){
    print("Please provide a file name")

  }else{
    saveRDS(object, file = paste0(path, "/", file, ".rds"))
  }
}

#' Save plot
#'
#' Wrapper function for saving plots as pdf.
#'
#' @param plot The plot.
#' @param path Directory path.
#' @param file The file name.
#' @param dim_h Height parameter for \code{\link[grDevices]{pdf}}
#' @param dim_w Width parameter for \code{\link[grDevices]{pdf}}

#' @return None
#'
#' @export
save.plot <- function(plot,
                      path=output_plot_path,
                      file,
                      dim_h=7,
                      dim_w=5){
  if(is.null(file)){
    print("Please provide a file name")

  }else{
    pdf(file = paste0(path, "/", file, ".pdf"), width = dim_w, height = dim_h)
    plot(plot)
    dev.off()
  }
}
