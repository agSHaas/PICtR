#' Data Import Helper 
#'
#' @param path The path to the directory containing the transformed and gated cytometry data files (string).  
#' @param filetype The file type of the cytometry data files (string, default: "csv").
#' @param sep The field separator character used in the cytometry data files. See \code{\link[data.table]{fread}} (string, default: `,`). 
#' @param verbose Verbosity (Boolean, default: `TRUE`).
#' @param downsample Optional random down-sampling of each individual file before they are combined into one data frame. If provided, the proportion of events that should be selected (numerical, default: `NULL`).
#'
#' @returns A `data.table`.
#' 
#' @export
import_cytometry_data <- function(path, 
                                  filetype = "csv", 
                                  sep = ",",
                                  verbose = T,
                                  downsample = NULL) {
  # find files 
  files <- dir(path, pattern=paste0("*\\.", filetype))
  
  # read values and save in list
  n <- length(files)
  file_list <- lapply(files, function(csv) {
    # read
    tab <- data.table::fread(paste0(path, csv), header = T, sep = sep)
    # extract meta data from file name 
    tab$sample <- csv
    
    if(!is.null(downsample)){
      tab <- dplyr::slice_sample(tab, prop = downsample)
    }
    
    if(verbose){
      message(paste("import", which(files == csv), "of", n, "done!\n")) 
    }
    return(tab)
  })
  
  # combine into one dataframe
  data <- dplyr::bind_rows(file_list)
  
  return(data)
}