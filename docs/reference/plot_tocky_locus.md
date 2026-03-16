# Generate basic QC plots for Tocky data (Timer-Blue vs Timer-Red 2d plots)

This function creates quick control plots for the TockyPrepData object
analyzing fluorescence changes over time in cellular activities.

## Usage

``` r
plot_tocky_locus(
  x,
  file = "PlotTockyLocus",
  n = 3,
  max_cell_number = 20000,
  viridis = FALSE,
  interactive = FALSE
)
```

## Arguments

- x:

  A TockyPrepData object produced by the function `prep_tocky`.

- file:

  The name of the output file.

- n:

  The number of plots per row and column in the output grid.

- max_cell_number:

  The maximum number of cells to be displayed per panel.

- viridis:

  (Optional). If TRUE, a colour-blind friendly colour set is used.

- interactive:

  (Optional). If TRUE, an interactive session is used to trim plot area.

## Value

An unchanged TockyPrepData object, primarily for consistency in pipeline
usage.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_tocky_locus(data)
} # }
```
