# Convert Timer Angle Data into color code

This function assigns colors to different ranges of angle values, with
an option to use colorblind-friendly colors from the viridis palette.

## Usage

``` r
Locus_to_colour(x, viridis = FALSE)
```

## Arguments

- x:

  Angle numeric vector.

- viridis:

  Logical, whether to use the viridis color palette.

## Value

a character vector for color code.

## Examples

``` r
if (FALSE) { # \dontrun{
col <- Locus_to_colour(x = c(0, 25, 45, 65, 90), viridis = TRUE)
} # }
```
