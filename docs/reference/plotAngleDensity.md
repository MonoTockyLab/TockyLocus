# Plot Density of Angles by Group Using Ridge Plots

This function takes a TockyPrepData object, which should have been
previously processed using the `timer_transform` function, and creates a
ridge plot showing the density distribution of angles for each group
defined in the dataset.

## Usage

``` r
plotAngleDensity(x, alpha = 0.3, group_order = NULL, scale = 2, legend = FALSE)
```

## Arguments

- x:

  A TockyPrepData object that has been processed with the
  `timer_transform` function.

- alpha:

  A number between 0 and 1 to be usedby ggridges.

- group_order:

  Optional. A character vector to define the order of group

- scale:

  A scaling factor to scale the height of the ridgelines. Used by
  ggridges.

- legend:

  Logical. If TRUE, legend is included.

## Value

A ggplot object showing the density distribution of angles by group.

## Examples

``` r
if (FALSE) { # \dontrun{
plotAngleDensity(x)
} # }
```
