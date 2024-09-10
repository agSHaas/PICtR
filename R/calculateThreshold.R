#' Threshold calculation using Otsu's method
#'
#' Calculates a threshold based on Ostu's method for a given histogram.
#'
#' @param hist A histogram to calculate the Otsu threshold for, for example created by \code{\link[graphics]{hist}}.
#' Please make sure that the histogram comprises enough breaks (n>100).
#'
#' @examples data("demo_lcmv")
#' demo_lcmv$ratio <- demo_lcmv$FSC.A/demo_lcmv$FSC.H
#' threshold <- calculateThreshold(hist(demo_lcmv$ratio, breaks=200))
#'
#' @return Numerical threshold value
#'
#'
#' @export
calculateThreshold <- function(hist){
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
