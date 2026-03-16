# Produce scatter plots of percentages of cells in each Tocky Locus.

Produce scatter plots of percentages of cells in each Tocky Locus.

## Usage

``` r
plotTockyLocus(
  x,
  percentTimer = FALSE,
  group_order = NULL,
  locus_colours = NULL,
  group_colors = NULL,
  group_by = TRUE,
  p_adjust_method = "fdr",
  ylim = NULL,
  stats = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  A TockyPrepData object

- percentTimer:

  A logical value for whether Percent Timer data is produced. Default is
  FALSE and produces Percent Parent data.

- group_order:

  The order of groups (optional).

- locus_colours:

  (optional) to choose colours for Tocky Loci.

- group_colors:

  (optional) to choose colours for groups.

- group_by:

  A logical value for whether different groups are plotted in different
  panels.

- p_adjust_method:

  A method for p-value adjustment in statistical tests.

- ylim:

  (Optional) the range of y values to be displayed.

- stats:

  A logical value for whether to produce statistical outputs. This is
  effective only for two-group analysis.

- verbose:

  Logical indicating whether to print Tocky Locus stats. Default is
  `TRUE`.

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
plotTockyLocus(x)
} # }
```
