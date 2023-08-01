save.obj <- function(object, path=output_obj_path, file){
  if(is.null(file)){
    print("Please provide a file name")
  }else{
    saveRDS(object, file = paste0(path, "/", file, ".rds"))
  }
}

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

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_facs <- function(object, marker1, marker2, var){
  object <- object
  index1 <- which(colnames(object)==marker1)
  index2 <- which(colnames(object)==marker2)
  object$density <- get_density(object[,index1], object[,index2], n = 100)
  p <- ggplot(object, aes(x=object[,index1], y=object[,index2], color=object[,var]))+ #color=density
    geom_point(size=1)+
    scale_color_gradientn(colours =pals::parula(1000))+
    theme_classic()+
    xlab(colnames(object)[index1])+
    ylab(colnames(object)[index2])+
    theme(axis.text = element_text(size=16),
          axis.title = element_text(size=16))
  return(p)
}

calculateThreshold<- function(hist){
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



sketch_wrapper <- function(channel=channel,
                          meta_data=NULL,
                          assay="FACS",
                          FSC.A="FSC.A.x", 
                          FSC.H="FSC.H.y",
                          n_sketch_cells=50000,
                          n_dims=dim(channel)[2], 
                          resolution=c(0.5,1,2,3,4),
                          obj_name="obj_sketched_non_projected",
                          group_by=NULL,
                          verbose=TRUE,
                          BPcell_dir=NULL,
                          working_dir=getwd()){
  if(packageVersion("Seurat") < "4.9.9.9058"){
    stop("Please install Seurat V5. If you already have Seurat V5 please install the most recent version from git")
  }else{
    options(future.globals.maxSize = 1e9)
    options(Seurat.object.assay.version = "v5")
  }
  
  if(dir.exists(paste0(working_dir, "/counts")) & is.null(BPcell_dir)){
    stop("Seems like BPcells was already run, either remove the folder 
            containing compressed data or define a path (BPcell_dir) to read in compressed files")
  }else if(dir.exists(paste0(working_dir, "/counts")) | !is.null(BPcell_dir)){
    warning("The pre-exisiting compressed data generated with BPCells will be used")
  }
  
  pkg <- c("Seurat", "tidyr","tidyverse", "dplyr", "BPCells", "readr")
  invisible(lapply(pkg, library, character.only = TRUE))
  
  obj <- CreateSeuratObject(as(object=t(channel), Class="dgCMatrix"), "FACS")
  
  if(dir.exists(paste0(working_dir, "/counts")) & is.null(BPcell_dir)){
    counts.mat <- open_matrix_dir(dir =  paste0(working_dir, "/counts"))
  }else if(!is.null(BPcell_dir)){
    counts.mat <- open_matrix_dir(dir =  BPcell_dir)
  }else{
    write_matrix_dir(mat = obj[["FACS"]]$counts, 
                     dir = paste0(working_dir, "/counts"))
    counts.mat <- open_matrix_dir(dir =  paste0(working_dir, "/counts"))
  }
  obj[["FACS"]]$counts <- counts.mat
  
  if(!is.null(meta_data)){
    obj@meta.data <- cbind(obj@meta.data, meta_data)
  }else{
    message("Please not that there are not metadata are added to the object")
  }
  obj$ratio <- obj@meta.data[,FSC.A]/obj@meta.data[,FSC.H]
  
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
    obj@meta.data <- left_join(rownames_to_column(obj@meta.data), 
                               dplyr::select(rownames_to_column(ratio.df), rowname, ratio_anno_group), by=c("rowname"="rowname")) %>% 
      column_to_rownames()
    obj$ratio_anno_group <- factor(obj$ratio_anno_group, levels = c("Ratio_low", "Ratio_high"))
   
  }
  obj@assays$FACS$data <- obj@assays$FACS$counts
  obj <- FindVariableFeatures(obj, verbose=verbose)
  obj <- SketchData(
    object = obj,
    ncells = n_sketch_cells,
    method = "LeverageScore",
    sketched.assay = "sketch", 
    verbose=verbose
  )
  
  DefaultAssay(obj) <- "sketch"
  obj <- obj %>% FindVariableFeatures(verbose=verbose) %>% 
    ScaleData(verbose=verbose) %>% 
    RunPCA(npcs=n_dims-1, verbose=verbose) %>%
    FindNeighbors(dims = 1:n_dims-1, verbose=verbose) %>%
    FindClusters(resolution = resolution, verbose = verbose) %>% 
    RunUMAP(features = 1:n_dims-1, return.model = TRUE, verbose = verbose)

  
  SeuratObject::saveRDS(object = obj,
                        file = paste0(obj_name, ".rds"),
                        destdir = working_dir)
  #unlink(paste0(working_dir, "counts"), recursive = TRUE)
  
  return(obj)
}

wrapper_for_plots <- function(obj=obj, 
                              feature_plot=TRUE, 
                              cluster_plot=TRUE, 
                              meta_list=list("ratio_anno"),
                              feature_plot_colors=pals::parula(1000),
                              ratio_plot_color=c(Ratio_low="dodgerblue2", Ratio_high="gold2"),
                              reduction="umap", 
                              label_size=3,
                              label_box=FALSE,
                              assay="sketch"){
  
  list.of.packages <-  c("cowplot", "pals", "khroma")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  pkg <- c("cowplot", "pals", "khroma")
  invisible(lapply(pkg, library, character.only = TRUE))
  smooth_rainbow <- colour("smooth rainbow")
  DefaultAssay(obj) <- assay
  # Feature Plots
  if(feature_plot){
    Feature_Plot <- FeaturePlot(obj, features=rownames(obj@assays$FACS), combine=FALSE, raster =TRUE, reduction = reduction)
    
    for(i in 1:length(Feature_Plot)) suppressMessages({
      Feature_Plot[[i]] <- Feature_Plot[[i]] + 
        scale_colour_gradientn(colours=feature_plot_colors) +
        #ggtitle(plot_names[i])+
        theme_classic()+
        NoLegend()+NoAxes()
    })
    f_nrow <- round(length(Feature_Plot)/8)
    
  }else{
    feature_plot <- "not caluculated"
  }
  # Cluster Plots
  if(cluster_plot){
    cluster_plots <- lapply(seq(1:length(grep("sketch_snn_res", colnames(obj@meta.data)))), function(i){
      name <- grep("sketch_snn_res", colnames(obj@meta.data), value = T)[i]
      cluster_plots <- DimPlot(object = obj,
                               group.by = name,
                               cols = smooth_rainbow(max(discard(as.numeric(as.character(obj@meta.data[,name])), is.na))+1, 
                                                     range = c(0.01, 0.99)), 
                               label = TRUE, label.box = label_box, label.size = label_size, repel = FALSE, reduction = reduction)
    })
    c_nrow <- round(length(grep("sketch_snn_res", colnames(obj@meta.data)))/3)
  }else{
    cluster_plot <- "not caluculated"
  }
  # Ratio_Plot
  if(any(meta_list=="ratio_anno")){
    ratio_plots <- DimPlot(object = obj,
                           group.by = "ratio_anno",
                           cols = ratio_plot_color, 
                           label = TRUE, repel = TRUE, label.size = label_size, reduction = reduction)
  }
  
  if(length(meta_list) > 0){
    meta_plots <-lapply(seq(1:length(meta_list)), function(i){
      name <- meta_list[[i]]
      if(is.numeric(obj@meta.data[,meta_list[[i]]])){
        FeaturePlot(obj, features=name, combine=T, raster =TRUE, reduction = reduction)+scale_colour_gradientn(colours=feature_plot_colors) 
      }else{
        if(length(unique(obj@meta.data[,meta_list[[i]]]))<=12){
          name <- meta_list[[i]]
          plot <- DimPlot(object = obj, group.by = name,
                          cols = pals::tol(12), 
                          label = TRUE, label.box = label_box,
                          label.size = label_size, repel = FALSE, 
                          reduction = reduction, combine = F)
        }else if(length(unique(obj@meta.data[,meta_list[[2]]]))<=25 && length(unique(obj@meta.data[,meta_list[[2]]]))>12){
          name <- meta_list[[i]]
          plot <- DimPlot(object = obj, group.by = name,
                          cols = pals::tol.rainbow(25), 
                          label = TRUE, label.box = label_box, 
                          label.size = label_size, repel = FALSE, 
                          reduction = reduction, combine = F)
        }else{
          name <- meta_list[[i]]
          cluster_plots <- DimPlot(object = obj,
                                   group.by = name,
                                   cols = smooth_rainbow(max(discard(as.numeric(as.character(obj@meta.data[,name])), is.na))+1, 
                                                         range = c(0.01, 0.99)), 
                                   label = TRUE, label.box = label_box,
                                   label.size = label_size, repel = FALSE,
                                   reduction = reduction, combine = F)
          
        }
      }
    })
    m_nrow <-  round((length(meta_list)/3))
    names(meta_plots) <- meta_list
  }else{
    meta_plots <- "not caluculated"
  }
  return(list(feature_plots= cowplot::plot_grid(plotlist = Feature_Plot, nrow = f_nrow),
              cluster_plots=cowplot::plot_grid(plotlist = cluster_plots, nrow = c_nrow),
              ratio_plot=ratio_plots, 
              meta_plots_grid=cowplot::plot_grid(plotlist = meta_plots, nrow = m_nrow),
              meta_plot_list=meta_plots))
}
