% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TockyLocus.R
\name{angle_to_locus}
\alias{angle_to_locus}
\title{Classify Angle Data into Tocky Loci}
\usage{
angle_to_locus(angle)
}
\arguments{
\item{angle}{A \verb{A numeric vector} of angles to classify.}
}
\value{
Tocky Locus data as a character vector.
}
\description{
This function classifies each value in a vector of angles into specific loci categories
based on the angle's value. It categorizes angles into 'New', 'NPt', 'Persistent',
'PAt', 'Arrested', or \code{NA} for Timer negative cells.
}
\examples{
\dontrun{
angle <- c(0, 29, 45, 65, 90)
locus <- angle_to_locus(angle)
}
}
\keyword{internal}
