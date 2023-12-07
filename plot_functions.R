save.obj <- function(object, path=output_obj_path, file){
  if(is.null(file)){
    print("Please provide a file name")
  
  }else{
    if (packageVersion("Seurat") >= "5.0.0"){
      SeuratObject::SaveSeuratRds(object = object,
                                  file = paste0(path, "/", file, ".rds"))
    }else{
      saveRDS(object, file = paste0(path, "/", file, ".rds"))
    }
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

wrapper_for_plots <- function(obj=obj, 
                              feature_plot=TRUE, 
                              cluster_plot=TRUE, 
                              meta_list=list("ratio_anno"),
                              feature_plot_colors=pals::parula(1000),
                              ratio_plot_color=c(Ratio_low="dodgerblue2", Ratio_high="gold2"),
                              reduction="umap", 
                              alpha=1,
                              raster=TRUE,
                              label_size=3,
                              cluster_handel="sketch_snn_res",
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
    Feature_Plot <- FeaturePlot(obj, features=rownames(obj[[assay]]), alpha = alpha, combine=FALSE, raster =raster, reduction = reduction)
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
    cluster_plots <- lapply(seq(1:length(grep(paste("^",cluster_handel, sep=""), colnames(obj@meta.data)))), function(i){
      name <- grep(cluster_handel, colnames(obj@meta.data), value = T)[i]
      cluster_plots <- DimPlot(object = obj,
                               group.by = name,
                               cols = smooth_rainbow(max(as.numeric(as.character(obj@meta.data[,name])), na.rm = T)+1, 
                                                     range = c(0.01, 0.99)), alpha = alpha, raster=raster,
                               label = TRUE, label.box = label_box, label.size = label_size, repel = FALSE, reduction = reduction)
    })
    c_nrow <- ifelse(ceiling((length(meta_list)/2))==0,1,ceiling((length(meta_list)/2)))
  }else{
    cluster_plot <- "not caluculated"
  }
  # Ratio_Plot
  if(any(meta_list=="ratio_anno")){
    ratio_plots <- DimPlot(object = obj,
                           group.by = "ratio_anno",
                           cols = ratio_plot_color,
                           alpha = alpha, raster=raster,
                           label = TRUE, label.box = label_box, label.size = label_size, repel = FALSE, reduction = reduction)
  }else{
    ratio_plots <- NULL
  }
  
  if(length(meta_list) > 0){
    meta_plots <-lapply(seq(1:length(meta_list)), function(i){
      name <- meta_list[[i]]
      #print(name)
      if(is.numeric(obj@meta.data[,meta_list[[i]]])){
        FeaturePlot(obj, features=name, combine=T, alpha = alpha, raster =TRUE, reduction = reduction)+scale_colour_gradientn(colours=feature_plot_colors) 
      }else{
        if(length(unique(obj@meta.data[,meta_list[[i]]]))<=12){
          name <- meta_list[[i]]
          plot <- DimPlot(object = obj, group.by = name,
                          cols = pals::tol(12), 
                          alpha = alpha, raster=raster,
                          label = TRUE, label.box = label_box,
                          label.size = label_size, repel = FALSE, 
                          reduction = reduction, combine = F)
        }else if(length(unique(obj@meta.data[,meta_list[[2]]]))<=25 && length(unique(obj@meta.data[,meta_list[[2]]]))>12){
          name <- meta_list[[i]]
          plot <- DimPlot(object = obj, group.by = name,
                          cols = pals::tol.rainbow(25), 
                          alpha = alpha, raster=raster,
                          label = TRUE, label.box = label_box, 
                          label.size = label_size, repel = FALSE, 
                          reduction = reduction, combine = F)
        }else if(is.factor(obj@meta.data[,name])){
          name <- meta_list[[i]]
          plot <- DimPlot(object = obj,
                          group.by = name,
                          cols = smooth_rainbow(max(discard(as.numeric(as.character(obj@meta.data[,name])), is.na))+1, 
                                                range = c(0.01, 0.99)), 
                          label = TRUE, label.box = label_box,
                          label.size = label_size, repel = FALSE,
                          alpha = alpha, raster=raster,
                          reduction = reduction, combine = F)
          
        }else{
          index <- length(unique(obj@meta.data[,name]))
          plot <- DimPlot(object = obj,
                          group.by = name,
                          cols = smooth_rainbow(index, range = c(0.01, 0.99)), 
                          label = TRUE, label.box = label_box,
                          label.size = label_size, repel = FALSE,
                          alpha = alpha, raster=raster,
                          reduction = reduction, combine = F)
        }
      }
    })
    m_nrow <-  ifelse(round((length(meta_list)/3))==0,1,round((length(meta_list)/3)))
    names(meta_plots) <- meta_list
  }else{
    meta_plots <- "not caluculated"
  }
  return(list(feature_plots= cowplot::plot_grid(plotlist = Feature_Plot, nrow = f_nrow),
              cluster_plots=cowplot::plot_grid(plotlist = cluster_plots, nrow = c_nrow),
              ratio_plots=ratio_plots, 
              meta_plots_grid=cowplot::plot_grid(plotlist = meta_plots, nrow = m_nrow),
              meta_plot_list=meta_plots))
}

split_plot_sketch <- function(obj, 
                              group_by="seurat_clusters", 
                              split_by="ratio_anno"){
  install.packages(setdiff("ggrastr", rownames(installed.packages())))
  cells <- rownames(obj@meta.data)[!is.na(obj$seurat_clusters)]
  meta <- obj@meta.data[cells,]
  meta <- cbind(meta, obj@reductions$umap@cell.embeddings[cells,])
  
  pkg <- c("cowplot", "pals", "khroma")
  invisible(lapply(pkg, library, character.only = TRUE))
  smooth_rainbow <- colour("smooth rainbow")
  p <- meta %>% ggplot(aes(umap_1, umap_2, color=meta[[group_by]]))+ggrastr::geom_point_rast(aes(alpha=0.1), size=0.3)+
    facet_wrap(~meta[[split_by]])+
    theme_classic()+
    guides(alpha = "none")+
    labs(color=group_by)+
    {if(length(unique(meta[[group_by]]))<=12){
      scale_color_manual(values = pals::tol(12))}
      else{
        scale_color_manual(values = smooth_rainbow(max(discard(as.numeric(as.character(obj@meta.data[,group_by])), is.na))+1, 
                                                   range = c(0.01, 0.99)))
      }}
  return(p)
}

ratio_cluster_plot <- function(obj, 
                               clusters="seurat_clusters", 
                               ratio="ratio_anno",
                               assay="FACS"){
  # load dependencies
  pkg <- c("ggplot2")
  invisible(lapply(pkg, library, character.only = TRUE))
  
  # switch to specified assay
  DefaultAssay(obj) <- assay
  # number of clusters 
  max_i <- max(as.numeric(as.vector(obj@meta.data[,clusters])), na.rm = T)
  original_warning <- options(warn = -1)
  # plot
  obj@meta.data %>%  group_by(.data[[clusters]], .data[[ratio]]) %>% count(.data[[ratio]]) %>% 
    ggplot(aes(x = reorder(.data[[clusters]], as.numeric(.data[[clusters]]), FUN = max), y = n, fill = .data[[ratio]])) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = rev(c("orange2", "dodgerblue2"))) +
    scale_x_discrete(limits = as.character(0:max_i)) +  # Reorder levels within the range of 1 to i
    xlab("Clusters") +
    labs(y = "Count") +
    ggtitle("Cluster Analysis") +
    theme_minimal()
}


# Code is taken from here:https://psyteachr.github.io/msc-conv/visualisation.html
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}


ClusterMarker <- function(obj, 
                          assay="FACS",
                          cluster_handel="cluster_full", 
                          #resolution=0.5,
                          sort=TRUE){
  #cluster_res <- paste0(cluster_handel, ".", resolution)
  dat <- setDT(as.data.frame(BPCells::as.matrix(obj[[assay]]$counts)), keep.rownames = "rowname")[] %>% 
    pivot_longer(-rowname) %>% 
    left_join(obj@meta.data %>% rownames_to_column() %>% dplyr::select(rowname, .data[[cluster_handel]]), by=c("name"="rowname"))
  
  p_list <- lapply(seq(1: max(as.numeric(as.vector(obj@meta.data[,cluster_handel]))+1)), function(i){
    
    dat$group <- ifelse(dat[,cluster_handel]==i-1, "case", "control")
    dat$col <- "dummy"
    
    dat <- dat %>% 
      group_by(rowname) %>% 
      mutate(mean_control=mean(value)) %>% 
      ungroup() %>% 
      group_by(rowname, .data[[cluster_handel]]) %>% 
      mutate(mean_case=mean(value)) 
    
    dat$col[dat$group=="control"]="gray80"
    dat$col[dat[,cluster_handel]==i-1 & dat$group=="case" & dat$mean_case >= dat$mean_control]="magenta4"
    dat$col[dat[,cluster_handel]==i-1 & dat$group=="case" & dat$mean_case <= dat$mean_control]="steelblue3"
    
    dat$rowname <- as.factor(dat$rowname)
    if(sort==TRUE){
      dat$diff <- dat$mean_control - dat$mean_case
      dat_new <- dat %>% dplyr::filter(.data[[cluster_handel]]==i-1) %>% 
        ungroup() %>% select(rowname, diff) %>% unique() %>% 
        arrange(diff)
      
      newFactor <- factor(dat$rowname, levels = dat_new$rowname)
      dat$rowname<-  newFactor 
    }
    
    p <- ggplot(dat, aes(x=rowname, y=value, fill= col))+
      geom_split_violin(trim = T, alpha = .4, scale = "width")+
      geom_boxplot(width = .2, alpha = .8, outlier.shape  = NA,
                   position = position_dodge(.25))+
      scale_fill_identity(guide = 'legend', labels=c("control", "Clst_high", "Clst_low"))+
      xlab("")+ylab("Normalized and transformed FMI")+
      ggtitle(paste0("Cluster_", i-1))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      guides(fill = guide_legend(reverse = T, title = ""))
    
    return(p)
  })
}


umap_domi_style <- function(data=obj,
                            group.by="seurat_clusters",
                            raster.dpi=500,
                            label=TRUE,
                            cols=pals::tol.rainbow(70),
                            umap1="umap_1",
                            umap2="umap_2",
                            reduction="umap",
                            raster=TRUE){
  if(isS4(data)){
    data$umap_1 <- data@reductions[[reduction]]@cell.embeddings[,1]
    data$umap_2 <- data@reductions[[reduction]]@cell.embeddings[,2]
    if(lable==TRUE){
      cluster_centers <- data@meta.data %>%
        group_by(.data[[group.by]]) %>%
        summarise(center_x = mean(umap_1), center_y = mean(umap_2))
    }
    if(raster==TRUE){
      data@meta.data %>% 
        ggplot(aes(umap_1, umap_2, color=.data[[group.by]])) +
        ggrastr::geom_point_rast(size=0.2,color="black", raster.dpi=700)+
        ggrastr::geom_point_rast(aes(color=as.factor(.data[[group.by]])),size=0.05, raster.dpi=700)+
        scale_color_manual(values =cols)+theme_classic()+
        guides(colour = guide_legend(override.aes = list(size=5)))+
        geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    }else(
      data@meta.data %>% 
        ggplot(aes(umap_1, umap_2, color=.data[[group.by]])) +
        geom_point(size=0.2,color="black")+
        geom_point(aes(color=as.factor(.data[[group.by]])),size=0.05)+
        scale_color_manual(values =cols)+theme_classic()+
        guides(colour = guide_legend(override.aes = list(size=5)))+
        geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    )
  }else if(is.data.frame(data)){
    if(!any(colnames(data)==umap1 | colnames(data)==umap2)){
      stop("Plase insert a umap_1 and umap_2 embedding values")
    }
    if(raster==TRUE){
      data %>% 
        ggplot(aes(umap_1, umap_2, color=.data[[group.by]])) +
        ggrastr::geom_point_rast(size=0.2,color="black", raster.dpi=700)+
        ggrastr::geom_point_rast(aes(color=as.factor(.data[[group.by]])),size=0.05, raster.dpi=700)+
        scale_color_manual(values =cols)+theme_classic()+
        guides(colour = guide_legend(override.aes = list(size=5)))+
        geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    }else(
      data %>% 
        ggplot(aes(umap_1, umap_2, color=.data[[group.by]])) +
        geom_point(size=0.2,color="black")+
        geom_point(aes(color=as.factor(.data[[group.by]])),size=0.05)+
        scale_color_manual(values =cols)+theme_classic()+
        guides(colour = guide_legend(override.aes = list(size=5)))+
        geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    )
  }
}

