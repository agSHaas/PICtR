#' Demo data set
#'
#' Demo data set for the PICtR workflow containing sampled cells from the spleens of LCMV infected mice (day 7) and control.
#' For demonstration purposes, only CD3, MHCII, CD11c, CD11b, CD45_2, CD19, Ly6G, CD90_1, CD4, and CD8 are included as cell type markers.
#'
#' @docType data
#'
#' @usage data(demo_lcmv)
#'
#' @format 72620 observations (cells) and 14 variables (markers, scaled to 0-1023).
#' \describe{
#'   \item{FSC.A}{Forward scatter area}
#'   \item{FSC.H}{Forward scatter height}
#'   \item{SSC.A}{Side scatter area}
#'   \item{SSC.H}{Side scatter height}
#'   \item{fluorescence markers}{CD3, MHCII, CD11c, CD11b, CD45_2, CD19, Ly6G, CD90_1, CD4, CD8}
#' }
#' @references Vonficht, Jopp-Saile, Yousefian, Flore et al. Ultra-high scale cytometry-based cellular interaction mapping. Nature Methods (2025)
"demo_lcmv"
