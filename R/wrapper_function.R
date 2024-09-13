#' PICtR workflow wrapper
#'
#' Characterizes cells into ratio_high and ratio_low cells for physically interacting cell analysis,
#' samples a representative subset of cells using \code{\link[Seurat]{SketchData}} and runs the standard Seurat analysis workflow.
#'
#' @param channel A data frame with the dimensionality cells x flow cytometry parameters.
#' @param meta_data A data frame with meta_data for every cell (Default=NULL).
#' @param assay A character string with the name of the assay to be created (Default="FACS").
#' @param FSC.A The name (string) of the column containing the FSC.A scatter parameter (Default="FACS.A").
#' @param FSC.H The name (string) of the column containing the FSC.H scatter parameter (Default="FACS.H").
#' @param n_sketch_cells The number of cells to be subsampled by \code{\link[Seurat]{SketchData}} (Default=50000).
#' @param n_components The number of components computed and used during analysis. "all" or given as an integer (Default="all").
#' @param resolution A numerical vector with the desired resolutions for clustering (Default: c(0.5,1,2,3,4)). Only used for Louvain, SLM, and Leiden clustering.
#' @param clst_algorithm Algorithm used for clustering. For Seurat's \code{\link[Seurat]{FindClusters}}: 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' Alternatively, "hdbscan" for Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN) on the UMAP embedding, see \code{\link[dbscan]{hdbscan}}),
#' "flowMeans" for non-parametric clustering and segmented-regression-based change point detection, see \code{\link[flowMeans]{flowMeans}},
#' or "flowSOM" for self-organizing maps, see \code{\link[Spectre]{run.flowsom}}, (Default: Louvain, 1)
#' @param min_clst_size minimum cluster size for HDBSCAN. See details at \code{\link[dbscan]{hdbscan}} (Default: 100).
#' @param meta.k meta k for FlowSOM clustering, see \code{\link[Spectre]{run.flowsom}} (Default: "auto").
#' @param obj_name The name used for storing the Seurat object (character).
#' @param ratio Boolean variable. If TRUE the FSC ratio (FSC.A/FSC.H) will be calculated.
#' @param thresholding_method Method for thresholding. One of "otsu", "triangle", "kmeans", or any method from the autothresholdr package. For details see ?calculateThreshold() (Default = "otsu").
#' @param hist_breaks Number of histogram breaks for methods that rely on histograms. Should be >100 (Default: 2000).
#' @param group_by Optional grouping parameter to calculate the FSC ratio (FSC.A/FSC.H) thresholding method for.
#' @param verbose Verbosity (Boolean).
#' @param BPcell_dir Optional directory with the counts matrix for \code{\link[BPCells]{open_matrix_dir}}.
#' @param overwrite Overwrite existing BPCells directory? (Default: FALSE)
#' @param working_dir Directory path to be used as working directory (character string).
#' @param seed Seed for reproducibility (Default: 42).
#'
#' @examples obj <- sketch_wrapper(channel = demo_lcmv,
#'                       meta_data = demo_lcmv,
#'                       n_sketch_cells = 5000,
#'                       clst_algorithm = 1,
#'                       ratio = TRUE,
#'                       thresholding_method = "otsu")
#'
#' @return Seurat object.
#'
#' @references Hao et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology (2023). doi: \url{https://doi.org/10.1038/s41587-023-01767-y}.
#' @references Van Gassen S et al. (2015) FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. Cytom Part J Int Soc Anal Cytol 87: 636-645. \url{https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625}.
#' @references Aghaeepour, N., Nikolic, R., Hoos, H.H. and Brinkman, R.R. (2011), Rapid cell population identification in flow cytometry data. Cytometry, 79A: 6-13. \url{https://doi.org/10.1002/cyto.a.21007}.
#' @references Campello, R.J.G.B., Moulavi, D., Sander, J. (2013). Density-Based Clustering Based on Hierarchical Density Estimates. In: Pei, J., Tseng, V.S., Cao, L., Motoda, H., Xu, G. (eds) Advances in Knowledge Discovery and Data Mining. PAKDD 2013. Lecture Notes in Computer Science(), vol 7819. Springer, Berlin, Heidelberg. \url{https://doi.org/10.1007/978-3-642-37456-2_14}
#'
#' @import Seurat
#'
#'
#' @export
sketch_wrapper <- function(channel=channel,
                           meta_data=NULL,
                           assay="FACS",
                           FSC.A="FSC.A",
                           FSC.H="FSC.H",
                           n_sketch_cells=50000,
                           n_components="all",
                           clst_algorithm=1,
                           resolution=c(0.5,1,2,3,4),
                           min_clst_size=100,
                           meta.k="auto",
                           obj_name="obj_sketched_non_projected",
                           ratio=TRUE,
                           thresholding_method="otsu",
                           hist_breaks = 2000,
                           group_by=NULL,
                           verbose=TRUE,
                           BPcell_dir=NULL,
                           overwrite=F,
                           working_dir=getwd(),
                           seed = 42){
  # check Seurat version
  if(packageVersion("Seurat") < "4.9.9.9058"){
    stop("Please install Seurat V5. If you already have Seurat V5 please install the most recent version from git")
  }else{
    options(future.globals.maxSize = 1e9)
    options(Seurat.object.assay.version = "v5")
  }

  # check if BP cells was already run
  if(dir.exists(paste0(working_dir, "/counts")) & is.null(BPcell_dir) & overwrite == F){
    stop("Seems like BPcells was already run, either remove the folder containing compressed data, define a path (BPcell_dir) to read in compressed files, or specify overwrite = TRUE")
  }else if(dir.exists(paste0(working_dir, "/counts")) | !is.null(BPcell_dir)){
    message("The pre-exisiting compressed data generated with BPCells will be used")
  } else if(overwrite == T) {
    message("The BPCells directory will be overwritten")
  }

  # calculate FSC area to height ratio
  if(ratio){
    message("Calculation of Ratio")
    channel["ratio"] <- channel[FSC.A]/channel[FSC.H]
    channel$ratio <- scales::rescale(as.numeric(channel$ratio), to = c(0, 1023))
  }

  # create Seurat obj and save counts on-disk with BPcells
  if (!file.exists(paste0(working_dir, obj_name, ".rds"))) {
    obj <- CreateSeuratObject(as(object=t(channel), Class="dgCMatrix"), assay = assay)

    if (dir.exists(paste0(working_dir, "/counts")) & !is.null(BPcell_dir)){
      counts.mat <- BPCells::open_matrix_dir(dir =  paste0(working_dir, "/counts"))
    }else if(!is.null(BPcell_dir)){
      counts.mat <- BPCells::open_matrix_dir(dir =  BPcell_dir)
    }else{
      BPCells::write_matrix_dir(mat = obj[[assay]]$counts,
                       dir = paste0(working_dir, "/counts"), overwrite = overwrite)
      counts.mat <- BPCells::open_matrix_dir(dir =  paste0(working_dir, "/counts"))
    }
    obj[[assay]]$counts <- counts.mat

    # save object according to Seurat version
    message("Your newly generated object will be saved under: ", paste0(working_dir, "/", obj_name, ".rds"))
    saveRDS(obj, file=paste0(working_dir, "/", obj_name, ".rds"))

  # load pre-existing object if present
  }else{
    message("Your generated object saved under ", paste0(working_dir, obj_name, ".rds"), " will be used")
    obj <- readRDS(paste0(working_dir, obj_name, ".rds"))

    if (packageVersion("Seurat") >= "5.0.0"){
      # explicitly load on-disk layer
      obj[[assay]]$counts <- BPCells::open_matrix_dir(dir = paste0(working_dir, "/counts"))
    }
  }

  # add meta data and add features to meta data as well
  if(!is.null(meta_data)){
    meta_data <- meta_data %>% dplyr::select(-(intersect(colnames(channel), colnames(meta_data)))) # remove cols present in both meta_data and channel
    obj@meta.data <- cbind(obj@meta.data, meta_data, channel) # add provided meta_data and all channel values to meta.data
    obj$ratio <- channel$ratio
  }else{
    obj@meta.data <- cbind(obj@meta.data, channel)
  }

  # calculate threshold and divide cells accordingly
  if(ratio){
    cutoff <- calculateThreshold(data = obj$ratio, method = thresholding_method, breaks = hist_breaks, seed = seed)
    obj$ratio_anno <- ifelse(obj$ratio>=cutoff, "Ratio_high", "Ratio_low")
    obj$ratio_anno <- factor(obj$ratio_anno, levels = c("Ratio_low", "Ratio_high"))

    if(!is.null(group_by)){
      ratio_list <-  split(obj@meta.data, f=obj@meta.data[,group_by])
      ratio_anno <- lapply(ratio_list, function(list){
        cutoff <- calculateThreshold(data = list$ratio, method = thresholding_method, breaks = hist_breaks, seed = seed)
        list$ratio_anno_group <- ifelse(list$ratio>=cutoff, "Ratio_high", "Ratio_low")
        return(list)
      })

      names(ratio_anno) <- NULL
      ratio.df <- do.call(rbind, ratio_anno)


      meta <- dplyr::left_join(tibble::rownames_to_column(obj@meta.data),
                                        dplyr::select(rownames_to_column(ratio.df), rowname, ratio_anno_group),
                                        by=c("rowname"="rowname")) %>% tibble::column_to_rownames()
      obj@meta.data <- meta
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
    verbose=verbose,
    seed = seed
  )

  message("Sketching is done")

  # standard Seurat workflow
  DefaultAssay(obj) <- "sketch"

  # number of components to be used
  if (n_components == "all") {
    n_dims <- dim(channel)[2]
  } else if (all.equal(n_components, as.integer(n_components))) {
    n_dims <- n_components
  } else {
    stop("Please provide an integer for the number of components to be used, or specify all")
  }

  obj <- obj %>%
    FindVariableFeatures(verbose=verbose) %>%
    ScaleData(verbose=verbose) %>%
    RunPCA(npcs=n_dims, approx=F, verbose=verbose) %>%
    FindNeighbors(dims = 1:n_dims, verbose=verbose)

  # clustering and UMAP
    if (clst_algorithm %in% c(1,2,3)) {
      obj <- obj %>%
        FindClusters(resolution = resolution, algorithm=clst_algorithm, verbose = verbose) %>%
        RunUMAP(dims = 1:n_dims, return.model = TRUE, verbose = verbose, seed.use = seed)

    } else if (clst_algorithm == 4) { # leiden via leidenalg
      obj <- obj %>%
        FindClusters(resolution = resolution, algorithm=clst_algorithm, verbose = verbose, method = "igraph") %>%
        RunUMAP(dims = 1:n_dims, return.model = TRUE, verbose = verbose, seed.use = seed)

    } else if (clst_algorithm == "hdbscan") { # HDBSCAN on UMAP embedding

      if (!requireNamespace("dbscan", quietly = TRUE)) {
        stop("Package \"dbscan\" must be installed to use HDBSCAN.")
      }

      obj <- obj %>%
        RunUMAP(dims = 1:n_dims, return.model = TRUE, verbose = verbose, seed.use = seed)

      # extract UMAP embedding
      obj$umap_1 <- obj@reductions$umap@cell.embeddings[,1]
      obj$umap_2 <- obj@reductions$umap@cell.embeddings[,2]

      # HDBSCAN
      set.seed(seed)
      hdbscan <- dbscan::hdbscan(obj@meta.data %>% tibble::rownames_to_column() %>% dplyr::filter(rowname %in% Cells(obj)) %>% dplyr::select(umap_1, umap_2), minPts = min_clst_size)
      message("The clustering contains ", max(hdbscan$cluster), " cluster(s). Noise points are indicated as 0.")
      hdbscan <- as.data.frame(hdbscan$cluster, row.names = Cells(obj))
      meta <- left_join(obj@meta.data %>% tibble::rownames_to_column(), hdbscan %>% tibble::rownames_to_column(), by = "rowname") %>%
        tibble::column_to_rownames() %>%
        dplyr::rename(hdbscan_clusters = `hdbscan$cluster`) %>%
        dplyr::mutate(hdbscan_clusters = factor(hdbscan_clusters))
      obj@meta.data <- meta

    } else if (clst_algorithm == "flowMeans") { # flowMeans with euclidean distances to circumvent singularity issues with Mahalanobis

      if (!requireNamespace("flowMeans", quietly = TRUE)) {
        stop("Package \"flowMeans\" must be installed to use FlowMeans.")
      }

      set.seed(seed)
      flowmeans <- flowMeans::flowMeans(obj@meta.data %>% tibble::rownames_to_column() %>% dplyr::filter(rowname %in% Cells(obj)), varNames = Features(obj), Standardize = FALSE, Mahalanobis = FALSE, iter.max = 50)
      flowmeans <- as.data.frame(flowmeans@Label, row.names = Cells(obj))
      meta <- left_join(obj@meta.data %>% tibble::rownames_to_column(), flowmeans %>% tibble::rownames_to_column(), by = "rowname") %>%
        tibble::column_to_rownames() %>%
        dplyr::rename(flowMeans_clusters = `flowmeans@Label`) %>%
        dplyr::mutate(flowMeans_clusters = factor(flowMeans_clusters))
      obj@meta.data <- meta

      obj <- obj %>%
        RunUMAP(dims = 1:n_dims, return.model = TRUE, verbose = verbose, seed.use = seed)

    } else if (clst_algorithm == "flowSOM") { # FlowSOM via Spectre

      if (!requireNamespace("Spectre", quietly = TRUE)) {
        stop("Package \"Spectre\" >= v1.1.0 must be installed to use FlowSOM.")
      }

      meta <- Spectre::run.flowsom(obj@meta.data %>% tibble::rownames_to_column() %>% dplyr::filter(rowname %in% Cells(obj)), use.cols = Features(obj), meta.k = meta.k, clust.seed = seed, meta.seed = seed)
      meta <- meta %>%
        mutate_at(vars("FlowSOM_cluster", "FlowSOM_metacluster"), ~ factor(.))
      obj@meta.data <- meta

      obj <- obj %>%
        RunUMAP(dims = 1:n_dims, return.model = TRUE, verbose = verbose, seed.use = seed)

    } else {
      stop("Please indicate a clustering algorithm from this list: 1 (Louvain), 2 (Louvain algorithm with multilevel refinement), 3 (SLM algorithm), 4 (Leiden), hdbscan, flowMeans, flowSOM")
    }


  message("The object will be updated and saved")

  # save object
  saveRDS(obj, file=paste0(working_dir, "/", obj_name, ".rds"))

  return(obj)
}
