# Threshold calculation

Calculates a threshold based on different methods.

## Usage

``` r
calculateThreshold(data, method = "otsu", breaks = 2000, seed = 42)
```

## Arguments

- data:

  The data to calculate a threshold for. Must be non-negative integers
  for autothresholdr methods, see below.

- method:

  The method used for thresholding. One of: "otsu", "triangle", "kmeans"
  or methods from the autothresholdr package, see
  [`auto_thresh`](https://rorynolan.github.io/autothresholdr/reference/auto_thresh.html)
  ("IJDefault", "Huang", "Huang2", "Intermodes", "IsoData", "Li",
  "MaxEntropy", "Mean", "MinErrorI", "Minimum", "Moments", "Otsu",
  "Percentile", "RenyiEntropy", "Shanbhag", "Triangle" and "Yen"). For
  autothresholdr methods, the data will be forced to integers (Default:
  "otsu").

- breaks:

  Number of histogram breaks for methods that rely on histograms. Should
  be \>100 (Default: 2000).

- seed:

  Seed for reproducibility (Default: 42).

## Value

Numerical threshold value

## References

N. Otsu, "A Threshold Selection Method from Gray-Level Histograms," in
IEEE Transactions on Systems, Man, and Cybernetics, vol. 9, no. 1, pp.
62-66, Jan. 1979, doi:
[10.1109/TSMC.1979.4310076](https://agshaas.github.io/PICtR/reference/10.1109/TSMC.1979.4310076).

Zack GW, Rogers WE, Latt SA. Automatic measurement of sister chromatid
exchange frequency. Journal of Histochemistry & Cytochemistry.
1977;25(7):741-753. doi:
[10.1177/25.7.70454](https://agshaas.github.io/PICtR/reference/10.1177/25.7.70454).

Hartigan, J. A., & Wong, M. A. (1979). Algorithm AS 136: A K-Means
Clustering Algorithm. Journal of the Royal Statistical Society. Series C
(Applied Statistics), 28(1), 100–108. <https://doi.org/10.2307/2346830>.

G. Landini, D.A. Randell, S. Fouad, and A. Galton. Automatic
thresholding from the gradients of region boundaries. Journal of
Microscopy, 265(2), 185-195.

## Examples

``` r
data("demo_lcmv")
demo_lcmv$ratio <- demo_lcmv$FSC.A/demo_lcmv$FSC.H
threshold <- calculateThreshold(demo_lcmv$ratio, method = "otsu", breaks = 2000)
```
