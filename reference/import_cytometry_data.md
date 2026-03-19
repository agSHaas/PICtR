# Data Import Helper

Data Import Helper

## Usage

``` r
import_cytometry_data(
  path,
  filetype = "csv",
  sep = ",",
  verbose = T,
  downsample = NULL
)
```

## Arguments

- path:

  The path to the directory containing the transformed and gated
  cytometry data files (string).

- filetype:

  The file type of the cytometry data files (string, default: "csv").

- sep:

  The field separator character used in the cytometry data files. See
  [`fread`](https://rdrr.io/pkg/data.table/man/fread.html) (string,
  default: \`,\`).

- verbose:

  Verbosity (Boolean, default: \`TRUE\`).

- downsample:

  Optional random down-sampling of each individual file before they are
  combined into one data frame. If provided, the proportion of events
  that should be selected (numerical, default: \`NULL\`).

## Value

A \`data.table\`.
