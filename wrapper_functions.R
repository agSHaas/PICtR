calculateThreshold<- function(hist){
  # OTSU threshold
  hn = hist$counts
  wB = 0
  wF = 0
  mB =0
  mF =0
  total =sum(hist$counts)
  sumB =0
  sum = 0
  between=0
  maxi =0
  for (i in 1:(length(hist$breaks)-1)){
    sum = sum + ( i * hn[i])
  }
  for (i in 1:(length(hist$breaks)-1)) {
    wB = wB + hn[i]
    if(wB == 0){
      next
    }
    wF = total - wB
    if (wF == 0){
      break
    }
    sumB = sumB + ( i * hn[i])
    mB = sumB / wB;
    mF = (sum - sumB) / wF
    between = wB * wF * (mB - mF)^2;
    if ( between >= maxi ) {
      threshold1 = hist$breaks[i]
      if ( between > maxi ) {
        threshold2 = hist$breaks[i]
      }
      maxi = between          
    }
    
  }
  return ((threshold1 + threshold2 ) / 2)
}



#Sample per group
sketch_wrapper <- function(channel=channel,
                           meta_data=NULL,
                           assay="FACS",
                           FSC.A="FSC.A", 
                           FSC.H="FSC.H",
                           n_sketch_cells=50000,
                           n_dims, 
                           resolution=c(0.5,1,2,3,4),
                           obj_name="obj_sketched_non_projected",
                           group_by=NULL,
                           verbose=TRUE,
                           BPcell_dir=NULL,
                           ratio=TRUE,
                           working_dir=getwd()){
  # check Seurat version
  if(packageVersion("Seurat") < "4.9.9.9058"){
    stop("Please install Seurat V5. If you already have Seurat V5 please install the most recent version from git")
  }else{
    options(future.globals.maxSize = 1e9)
    options(Seurat.object.assay.version = "v5")
  }
  
  # check if BP cells was already run 
  if(dir.exists(paste0(working_dir, "/counts")) & is.null(BPcell_dir)){
    stop("Seems like BPcells was already run, either remove the folder containing compressed data or define a path (BPcell_dir) to read in compressed files")
  }else if(dir.exists(paste0(working_dir, "/counts")) | !is.null(BPcell_dir)){
    warning("The pre-exisiting compressed data generated with BPCells will be used")
  }
  
  # load dependencies
  pkg <- c("Seurat", "tidyr","tidyverse", "dplyr", "BPCells", "readr")
  invisible(lapply(pkg, library, character.only = TRUE))
  
  # calculate FSC area to height ratio 
  if(ratio){
    message("Calculation of Ratio")
    channel["ratio"] <- channel[FSC.A]/channel[FSC.H]
    #channel$ratio <- as.numeric(channel$ratio)
    channel$ratio <- scales::rescale(as.numeric(channel$ratio), to = c(0, 1023))
  }
  
  # create Seurat obj and save counts on-disk with BPcells 
  if(!file.exists(paste0(working_dir, obj_name, ".rds"))){
    obj <- CreateSeuratObject(as(object=t(channel), Class="dgCMatrix"), assay = assay)
    
    if(dir.exists(paste0(working_dir, "/counts")) & is.null(BPcell_dir)){
      counts.mat <- open_matrix_dir(dir =  paste0(working_dir, "counts"))
    }else if(!is.null(BPcell_dir)){
      counts.mat <- open_matrix_dir(dir =  BPcell_dir)
    }else{
      write_matrix_dir(mat = obj[[assay]]$counts, 
                       dir = paste0(working_dir, "/counts"), overwrite = TRUE)
      counts.mat <- open_matrix_dir(dir =  paste0(working_dir, "/counts"))
    }
    obj[[assay]]$counts <- counts.mat
    
    # save object according to Seurat version
    message("Your newly generated object will be saved under: ", paste0(working_dir, "/", obj_name, ".rds"))
    if (packageVersion("Seurat") >= "5.0.0"){
      SeuratObject::SaveSeuratRds(object = obj,
                                  file = paste0(working_dir, obj_name, ".rds"),
                                  destdir = paste0(working_dir, "counts/"))
    }else{
      SeuratObject::saveRDS(object = obj,
                            file = paste0(obj_name, ".rds"),
                            destdir = working_dir)
    }
  
  # load pre-existing object if present  
  }else{
    message("Your generated object saved under ", paste0(working_dir, obj_name, ".rds"), " will be used")
    obj <- readRDS(paste0(working_dir, obj_name, ".rds"))
    
    if (packageVersion("Seurat") >= "5.0.0"){
      # explicitly load on-disk layer
      obj[[assay]]$counts <- open_matrix_dir(dir = paste0(working_dir, "/counts"))
    }  
  }
  
  # add meta data
  if(!is.null(meta_data)){
    obj@meta.data <- cbind(obj@meta.data, meta_data)
    obj$ratio <- channel$ratio
  }else{
    message("Please note that there was no meta data added to the object")
  }
  
  # calculate threshold and divide cells accordingly 
  if(ratio){
    cutoff <- calculateThreshold(hist(obj$ratio, breaks = 2000, plot = FALSE))
    obj$ratio_anno <- ifelse(obj$ratio>=cutoff, "Ratio_high", "Ratio_low")
    obj$ratio_anno <- factor(obj$ratio_anno, levels = c("Ratio_low", "Ratio_high"))
    
    if(!is.null(group_by)){
      ratio_list <-  split(obj@meta.data, f=obj@meta.data[,group_by])  
      ratio_anno <- lapply(ratio_list, function(list){
        cutoff <- calculateThreshold(hist(list$ratio, breaks = 2000, plot = FALSE))
        list$ratio_anno_group <- ifelse(list$ratio>=cutoff, "Ratio_high", "Ratio_low")
        return(list)
      })
      
      names(ratio_anno) <- NULL
      ratio.df <- do.call(rbind, ratio_anno)
      
      
      obj@meta.data <- dplyr::left_join(tibble::rownames_to_column(obj@meta.data), 
                                        dplyr::select(rownames_to_column(ratio.df), rowname, ratio_anno_group), 
                                        by=c("rowname"="rowname")) %>% 
        tibble::column_to_rownames()
      obj$ratio_anno_group <- factor(obj$ratio_anno_group, levels = c("Ratio_low", "Ratio_high"))
    }
  }
  
  # sketching workflow
  message("Sketching started...")
  obj[[assay]]$data <- obj[[assay]]$counts
  obj <- FindVariableFeatures(obj, verbose=verbose)
  obj <- SketchData(
    object = obj,
    ncells = n_sketch_cells,
    method = "LeverageScore",
    sketched.assay = "sketch", 
    verbose=verbose
  )
  
  # standard Seurat workflow
  DefaultAssay(obj) <- "sketch"
  n_dims <- dim(channel)[2]
  obj <- obj %>% 
    FindVariableFeatures(verbose=verbose) %>% 
    ScaleData(verbose=verbose) %>% 
    RunPCA(npcs=n_dims-1, approx=F, verbose=verbose) %>%
    FindNeighbors(dims = 1:n_dims-1, verbose=verbose) %>%
    FindClusters(resolution = resolution, verbose = verbose) %>% 
    RunUMAP(dims = 1:n_dims-1, return.model = TRUE, verbose = verbose)
  
  message("Sketching is done")
  message("The object will be updated and saved")
  
  # save object according to Seurat version
  if (packageVersion("Seurat") >= "5.0.0"){
    SeuratObject::SaveSeuratRds(object = obj,
                                file = paste0(working_dir, obj_name, ".rds"),
                                destdir = paste0(working_dir, "counts/"))
  }else{
    SeuratObject::saveRDS(object = obj,
                          file = paste0(obj_name, ".rds"),
                          destdir = working_dir)
  }   
  #unlink(paste0(working_dir, "counts"), recursive = TRUE)
  return(obj)
}