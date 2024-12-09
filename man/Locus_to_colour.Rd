% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tocky_locus.R
\name{Locus_to_colour}
\alias{Locus_to_colour}
\title{Convert Timer Angle Data into color code}
\usage{
Locus_to_colour(x, viridis = FALSE)
}
\arguments{
\item{x}{Angle numeric vector.}

\item{viridis}{Logical, whether to use the viridis color palette.}
}
\value{
a character vector for color code.
}
\description{
This function assigns colors to different ranges of angle values, with an option
to use colorblind-friendly colors from the viridis palette.
}
\examples{
\dontrun{
col <- Locus_to_colour(x = c(0, 25, 45, 65, 90), viridis = TRUE)
}

}
