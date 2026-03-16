# Create cellcount table for calculating percentage parent

Create cellcount table for calculating percentage parent

## Usage

``` r
cellcount_table(x, percentTimer = TRUE)
```

## Arguments

- x:

  A data frame containing the Timer Angle variable with the column name
  Angle.

- percentTimer:

  A logical value for whether Percent Timer data is produced. Default is
  FALSE and produces Percent Parent data.

## Value

A data frame for the cell number of each Tocky Locus. NA is returned for
Timer negative cells.

## Examples

``` r
if (FALSE) { # \dontrun{
cellcounttable <- cellcount_table(x, percentTimer = FALSE)
} # }
```
