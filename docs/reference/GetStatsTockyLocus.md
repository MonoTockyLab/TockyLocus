# Perform Statistical Tests for Tocky Locus Analysis

This function performs statistical tests on Tocky Locus data, allowing
for different methods and p-value adjustments.

## Usage

``` r
GetStatsTockyLocus(
  x,
  percentTimer = FALSE,
  p_adjust_method = "BH",
  method = "ASR",
  verbose = TRUE
)
```

## Arguments

- x:

  A `TockyPrepData` object containing Tocky Locus data.

- percentTimer:

  Logical. If `TRUE`, the percentages of Timer-positive cells will be
  used; if `FALSE`, percentages of parent cells will be used.

- p_adjust_method:

  Character string specifying the method for p-value adjustment in
  multiple testing. Default is `'BH'` (Benjamini-Hochberg). Other
  methods available in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html), such as `'holm'`
  or `'bonferroni'`, can also be used.

- method:

  Character string specifying the statistical test method to use.
  Options are:

  `'Wilcox'`

  :   Mann-Whitney U test (Wilcoxon rank sum test) without data
      transformation.

  `'ASR'`

  :   Arcsine Square Root Transformation, followed by a normality test
      and t-test.

  `'Logit'`

  :   Logit Transformation, followed by a normality test and t-test.

- verbose:

  Logical indicating whether to print progress messages and outputs.
  Default is `TRUE`.

## Value

A `TockyPrepData` object containing the statistical outputs for Tocky
Locus Analysis, stored in `x@Tocky$TockyLocusStats`.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- GetStatsTockyLocus(x, method = 'ASR')
} # }
```
