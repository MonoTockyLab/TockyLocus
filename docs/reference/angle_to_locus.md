# Classify Angle Data into Tocky Loci

This function classifies each value in a vector of angles into specific
loci categories based on the angle's value. It categorizes angles into
'New', 'NPt', 'Persistent', 'PAt', 'Arrested', or `NA` for Timer
negative cells.

## Usage

``` r
angle_to_locus(angle)
```

## Arguments

- angle:

  A `A numeric vector` of angles to classify.

## Value

Tocky Locus data as a character vector.

## Examples

``` r
if (FALSE) { # \dontrun{
angle <- c(0, 29, 45, 65, 90)
locus <- angle_to_locus(angle)
} # }
```
