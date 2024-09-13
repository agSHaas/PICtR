#' Threshold calculation
#'
#' Calculates a threshold based on  method for a given histogram.
#'
#' @param data The data to calculate a threshold for. Must be non-negative integers for autothresholdr methods, see below.
#' @param method The method used for thresholding. One of: "otsu", "triangle", "kmeans" or methods from the autothresholdr package, see \code{\link[autothresholdr]{auto_thresh}} ("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean", "MinErrorI", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle" and "Yen"). For autothresholdr methods, the data will be forced to integers (Default: "otsu").
#' @param breaks Number of histogram breaks for methods that rely on histograms. Should be >100 (Default: 2000).
#' @param seed Seed for reproducibility (Default: 42).
#'
#'
#' @examples data("demo_lcmv")
#' demo_lcmv$ratio <- demo_lcmv$FSC.A/demo_lcmv$FSC.H
#' threshold <- calculateThreshold(demo_lcmv$ratio, method = "otsu", breaks = 2000)
#'
#' @return Numerical threshold value
#'
#' @references N. Otsu, "A Threshold Selection Method from Gray-Level Histograms," in IEEE Transactions on Systems, Man, and Cybernetics, vol. 9, no. 1, pp. 62-66, Jan. 1979, doi: \url{10.1109/TSMC.1979.4310076}.
#' @references Zack GW, Rogers WE, Latt SA. Automatic measurement of sister chromatid exchange frequency. Journal of Histochemistry & Cytochemistry. 1977;25(7):741-753. doi: \url{10.1177/25.7.70454}.
#' @references Hartigan, J. A., & Wong, M. A. (1979). Algorithm AS 136: A K-Means Clustering Algorithm. Journal of the Royal Statistical Society. Series C (Applied Statistics), 28(1), 100–108. \url{https://doi.org/10.2307/2346830}.
#' @references G. Landini, D.A. Randell, S. Fouad, and A. Galton. Automatic thresholding from the gradients of region boundaries. Journal of Microscopy, 265(2), 185-195.
#'
#'
#' @export
calculateThreshold <- function(data,
                               method = "otsu",
                               breaks = 2000,
                               seed = 42)
  {
  # OTSU THRESHOLD ----
  if (method == "otsu") {

    # for backwards compatibility
    if (class(data) == "histogram") {
      hist <- data
    } else {
      hist <- hist(data, breaks = breaks, plot = F)
    }

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

  # TRIANGLE ALGORITHM ----
  } else if (method == "triangle") {
    # code ported from ImageJ's implementation of the triangle algorithm in java
    # initialise variables
    min <- 0
    max <- 0
    min2 <- 0
    min3 <- 0
    dmax <- 0
    threshold <- 0
    hist <- hist(data, breaks = breaks, plot = F)

    # counts and breaks of the histogram
    counts <- hist$counts
    breaks <- hist$breaks

    # find first non-zero (minimum) bin of the histogram
    for (i in seq_along(counts)) {
      if (counts[i] > 0) {
        min <- i
        break
      }
    }
    # move one step back to the first zero point before min
    # note 1-based indexing in R
    if (min > 1) {
      min <- min - 1
    }

    # find last non-zero (minimum) bin of the histogram
    for (i in seq_along(counts)) {
      if (counts[length(counts) - i + 1] > 0) {
        min2 <- length(counts) - i + 1
        break
      }
    }
    # Move one step forward to the first zero point after min2
    # note 1-based indexing in R
    if (min2 < length(counts)) {
      min2 <- min2 + 1
    }


    # find the peak
    for (i in seq_along(counts)) {
      if (counts[i] > dmax) {
        max <- i
        dmax <- counts[i]
      }
    }

    # two possible thresholds for the two sides of the histogram - find the side
    # where the distance of the peak to the minumun is furthest
    inverted <- FALSE

    if ((max - min) < (min2 - max)) {
      # Reverse the histogram
      inverted <- TRUE
      counts <- rev(counts)
      min <- length(counts) - min2 + 1
      max <- length(counts) - max + 1
    }

    # If min and max are the same, return min
    if (min == max) {
      return(min)
    }

    # describe the line from the peak of the histogram to the minimum
    nx <- counts[max] # peak frequency
    ny <- min - max
    d <- sqrt(nx^2 + ny^2)
    nx <- nx / d
    ny <- ny / d
    d <- nx * min + ny * counts[min]

    # find the point with the maximum distance from the line connecting peak
    # and minimum to the histogram
    threshold <- min # initialize
    distance <- 0

    for (i in (min + 1):max) {
      new_distance <- nx * i + ny * counts[i] - d
      if (new_distance > distance) {
        threshold <- i
        distance <- new_distance
      }
    }

    # Adjust the split point
    threshold <- threshold - 1

    # reverse the histogram / threshold back if needed
    if (inverted) {
      threshold <- length(counts) - threshold + 1
    }

    return(breaks[threshold])

  # KMEANS CLUSTERING
  } else if (method == "kmeans") {
    set.seed(seed)

    km <- kmeans(data, centers = 2) # kmeans with 2 centers for binary classification
    return(mean(km$centers)) # threshold as the mean of the cluster centers

  # AUTOTHRESHOLDR METHODS
  } else if (method %in% c("IJDefault",
                           "Huang",
                           "Huang2",
                           "Intermodes",
                           "IsoData",
                           "Li",
                           "MaxEntropy",
                           "Mean",
                           "MinErrorI",
                           "Minimum",
                           "Moments",
                           "Otsu",
                           "Percentile",
                           "RenyiEntropy",
                           "Shanbhag",
                           "Triangle",
                           "Yen")) {

    if (!requireNamespace("autothresholdr", quietly = TRUE)) {
      stop("Package \"autothresholdr\" must be installed to use this thresholding method.")
    }

    threshold <- autothresholdr::auto_thresh(as.integer(data), method = method)[1]
    return(threshold)
  }
}
