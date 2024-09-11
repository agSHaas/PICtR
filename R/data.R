#' Demo data set
#'
#' Demo data set for the PICtR workflow containing 10 percent randomly sampled cells from the spleen of one LCMV infected mouse (day 7, Figure 5 and 6).
#' Data was analysed using high parametric spectral flow cytometry (n=26 markers).
#'
#' @docType data
#'
#' @usage data(demo_lcmv)
#'
#' @format 161287 observations (cells) and 30 variables (markers, scaled to 0-1023).
#' \describe{
#'   \item{FSC.A}{Forward scatter area}
#'   \item{FSC.H}{Forward scatter height}
#'   \item{SSC.A}{Side scatter area}
#'   \item{SSC.H}{Side scatter height}
#'   \item{fluorescence markers}{TCRb, CD3, NK11, MHCII, CD24, CD11c, ICOS, CD11b, CD451, CD19, CD117, Ly6G, CD138, Ly6C, CD901, CD25, CD452, CD68, F480, CD357, CD172, CD317, TCRgd, CD8, CD4, IgD}
#' }
#' @references Vonficht, Jopp-Saile, Yousefian, Flore et al. (2024). Ultra-high scale cytometry-based cellular interaction mapping.
"demo_lcmv"

#' Demo Seurat object
#'
#' Demo Seurat object based on the demo data set for the spleen of one LCMV infected mouse at day 7.
#'
#' @docType data
#'
#' @usage data(demo_obj)
#'
#' @format A Seurat object with 161287 samples (cells) and 31 Features. n = 5000 cells were sketched.
#' \describe{
#'   \item{meta.data}{Contains all scatter and fluorescence markers,
#'   the calculated FSC.A/FSC.H ratio,
#'   the classification into ratio_high/ratio_low cells (ratio_anno),
#'   the clustering solutions for sketched cells with resolutions of 0.8 and 1 (sketch_snn_res),
#'   and predicted cluster labels for all cells (clusters_predicted_obj).}
#'   \item{assay "FACS"}{The Seurat assay for n = 161287 cells.}
#'   \item{assay "sketch"}{The Seurat assay for n = 5000 sketched cells.}
#'   \item{reductions}{Calculated dimensional reductions (PCA, UMAP).}
#'   \item{misc}{List containing the selected clusters representing physically interacting cells.}
#' }
#'
#' @references Vonficht, Jopp-Saile, Yousefian, Flore et al. (2024). Ultra-high scale cytometry-based cellular interaction mapping.
"demo_obj"
