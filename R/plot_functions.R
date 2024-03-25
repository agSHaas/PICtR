#' Plot wrapper function.
#'
#' A wrapper function for plots that are often used during exploratory analysis within the \link[=https://satijalab.org/seurat/]{Seurat framework}.
#' Includes \code{\link[Seurat]{FeaturePlot}}, \code{\link[Seurat]{DimPlot}} for different parameters, and \code{\link{ratio_cluster_plot}}.
#'
#' @param obj The Seurat object.
#' @param feature_plot Boolean. TRUE indicates that UMAP is colored for all features (see \code{\link[Seurat]{FeaturePlot}})
#' @param cluster_plot Boolean. TRUE indicates that UMAP is colored by the different cluster resolutions stored with the cluster_handle.
#' @param meta_list List of meta.data columns to color \code{\link[Seurat]{DimPlot}} by.
#' Can be both numeric to generate \code{\link[Seurat]{FeaturePlot}} or class character or factor for \code{\link[Seurat]{DimPlot}}.
#' @param cluster_handle Prefix for the clustering solutions in the meta.data slot.
#' @param feature_plot_colors Color palette for \code{\link[Seurat]{FeaturePlot}}.
#' @param ratio_plot_color Color palette for \code{\link{ratio_cluster_plot}}.
#' @param reduction Reduction to use for plotting, for example UMAP.
#' @param alpha Alpha value for plotting.
#' @param raster Convert points to raster format. If TURE plot is rasterized to raster.dpi=c(512, 512).
#' @param label_size Size of the labels plotted within the embedding.
#' @param label_box Plot boxes around labels in the color of the cluster.
#' @param assay The Seurat assay (default FACS).
#'
#' @return A list with the requested plots.
#'
#' @importFrom pals parula
#' @importFrom khroma colour
#' @importFrom cowplot plot_grid
#' @importFrom grDevices dev.off pdf
#' @importFrom ggrastr geom_point_rast
#'
#' @export
wrapper_for_plots <- function(obj=obj,
                              feature_plot=TRUE,
                              cluster_plot=TRUE,
                              meta_list=list("ratio_anno"),
                              cluster_handle="sketch_snn_res",
                              feature_plot_colors=pals::parula(1000),
                              ratio_plot_color=c(Ratio_low="dodgerblue2", Ratio_high="gold2"),
                              reduction="umap",
                              alpha=1,
                              raster=TRUE,
                              label_size=3,
                              label_box=FALSE,
                              assay="sketch"){

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
    feature_plot <- "not calculated"
  }

  # Cluster Plots
  if(cluster_plot){
    cluster_plots <- lapply(seq(1:length(grep(paste("^",cluster_handle, sep=""), colnames(obj@meta.data)))), function(i){
      name <- grep(cluster_handle, colnames(obj@meta.data), value = T)[i]
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
    meta_plots <- "not calculated"
  }
  return(list(feature_plots= cowplot::plot_grid(plotlist = Feature_Plot, nrow = f_nrow),
              cluster_plots=cowplot::plot_grid(plotlist = cluster_plots, nrow = c_nrow),
              ratio_plots=ratio_plots,
              meta_plots_grid=cowplot::plot_grid(plotlist = meta_plots, nrow = m_nrow),
              meta_plot_list=meta_plots))
}

#' Split plot wrapper.
#'
#' Dimensional reduction (UMAP) plot split by a given parameter.
#' Per default the returned split plots are rasterized using \code{\link[ggrastr]{geom_point_rast}}.
#'
#' @param obj The Seurat object.
#' @param group_by Parameter to group the plot by.
#' @param split_by Parameter to split the plot by.
#'
#' @return ggplot object
#'
#' @importFrom pals parula
#' @importFrom khroma colour
#'
#' @export
split_plot_sketch <- function(obj,
                              group_by="seurat_clusters",
                              split_by="ratio_anno"){
  cells <- rownames(obj@meta.data)[!is.na(obj$seurat_clusters)]
  meta <- obj@meta.data[cells,]
  meta <- cbind(meta, obj@reductions$umap@cell.embeddings[cells,])

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

#' Ratio cluster plot.
#'
#' Stacked bar plot of each cluster with the proportion of cells below/above the threshold determined with Otsu's method using the FSC.A/FSC.H ratio.
#'
#' @param obj The Seurat object.
#' @param clusters The string of the meta.data column with the clustering resolution to plot.
#' @param ratio The meta.data column with the classification of cells (ratio_high/ratio_low) determined using the FSC.A/FSC.H ratio and the determined threshold using Otsu's method.
#' @param assay The Seurat assay to use (default FACS).
#'
#' @return None
#'
#' @importFrom stats reorder
#'
#' @export
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

#' Rasterized UMAP
#'
#' Wrapper function to plot a dimensional reduction plot colored by clusters.
#'
#' @param data Seurat object.
#' @param group.by meta.data column to group the dimensional reduction plot by, for example a clustering solution.
#' @param raster.dpi Pixel resolution (numeric; default 500).
#' @param label Boolean. Plot labels?
#' @param cols Color palette.
#' @param umap1 First UMAP dimension as found in the meta.data.
#' @param umap2 Second UMAP dimension as found in the meta.data.
#' @param reduction Which reduction to use.
#' @param raster Boolean. Rasterize the plot?
#'
#' @return Plot
#'
#' @importFrom ggrastr geom_point_rast
#'
#' @export
umap_rasterized <- function(data=obj,
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
    if(label==TRUE){
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
        ggrepel::geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    }else(
      data@meta.data %>%
        ggplot(aes(umap_1, umap_2, color=.data[[group.by]])) +
        geom_point(size=0.2,color="black")+
        geom_point(aes(color=as.factor(.data[[group.by]])),size=0.05)+
        scale_color_manual(values =cols)+theme_classic()+
        guides(colour = guide_legend(override.aes = list(size=5)))+
        ggrepel::geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
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
        ggrepel::geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    }else(
      data %>%
        ggplot(aes(umap_1, umap_2, color=.data[[group.by]])) +
        geom_point(size=0.2,color="black")+
        geom_point(aes(color=as.factor(.data[[group.by]])),size=0.05)+
        scale_color_manual(values =cols)+theme_classic()+
        guides(colour = guide_legend(override.aes = list(size=5)))+
        ggrepel::geom_text_repel(data = cluster_centers, aes(x = center_x, y = center_y, label = .data[[group.by]]), color="black", size = 3)
    )
  }
}

#' Marker Enrichment Modeling (MEM) Heatmap.
#'
#' @param obj The Seurat object.
#' @param markers Meta.data columns with features that should be plotted in the heat map and the clustering resolution.
#' @param cluster_col Character string specifying the column that contains the clustering solution.
#' @param cols Color palette.
#' @param heatmap_name Title of the heatmap.
#' @param heatmap_column_title Title for the columns of the heatmap.
#' @param heatmap_row_title Title for the rows of the heatmap.
#' @param scale_width Scaling factor for the width of the heatmap in relation to the number of columns.
#' @param scale_height Scaling factor for the height of the heatmap in relation to the number of rows.
#'
#' @return MEM heat map.
#'
#' @importFrom cytoMEM MEM
#' @importFrom pals coolwarm
#' @importFrom ComplexHeatmap Heatmap
#'
#' @export
MEM_heatmap <- function(obj,
                        markers = c(),
                        cluster_col = "seurat_clusters",
                        cols = pals::coolwarm(100),
                        heatmap_name = "MEM enrichment score",
                        heatmap_column_title = "marker",
                        heatmap_row_title = "cluster",
                        scale_width = 2.2,
                        scale_height = 5){

  # select markers and cluster annotations
  MEM <- obj@meta.data %>%
    dplyr::select(markers) %>%
    rename(cluster = .data[[cluster_col]]) %>%
    mutate(cluster = as.numeric(as.character(cluster)))

  # calculate MEM scores
  MEM_values <- cytoMEM::MEM(MEM,
                             transform = F,
                             choose.markers = F,
                             markers = "all",
                             choose.ref = F,
                             zero.ref = F,
                             IQR.thresh = NULL)

  # extract MEM matrix for plotting
  heatmap <- as.data.frame(MEM_values$MEM_matrix[[1]])

  # create heatmap
  h1 <- Heatmap(t(as.matrix(heatmap)),
                col=cols,
                name = "MEM enrichment score",
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D2",
                column_title = "marker",
                row_title = "cluster",
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 9),
                show_row_names = T,
                cluster_rows = T,
                cluster_columns = T,
                show_parent_dend_line = FALSE,
                width = ncol(heatmap)*unit(scale_width, "mm"),
                height = nrow(heatmap)*unit(scale_height, "mm"))

  # plot
  print(h1)
}
