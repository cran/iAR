% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_unidata.R
\name{unidata}
\alias{unidata}
\title{Unidata Class}
\usage{
unidata(
  times = integer(0),
  series = integer(0),
  series_esd = integer(0),
  series_names = character(0)
)
}
\arguments{
\item{times}{A numeric vector representing the time points.}

\item{series}{A numeric vector representing the values of the time series.}

\item{series_esd}{A numeric vector representing the error standard deviations of the time series.}

\item{series_names}{An optional character vector of length 1 representing the name of the series.}
}
\description{
The `unidata` class is an S7 class designed to represent univariate irregularly observed time series models
with associated times, values, and optional error standard deviations.
}
\section{Validation Rules}{

- `@times`, `@series`, and `@series_esd` must be numeric vectors.
- `@times` must not contain `NA` values and must be strictly increasing.
- The length of `@series` must match the length of `@times`.
- The length of `@series_esd` must be 0, 1, or equal to the length of `@series`.
- `NA` values in `@series` must correspond exactly (positionally) to `NA` values in `@series_esd`.
- `@series_names`, if provided, must be a character vector of length 1.
}

\examples{
# Create a unidata object
unidata_instance <- unidata(
  times = c(1, 2, 3, 4),
  series = c(10, 20, 15, 25),
  series_esd = c(1, 1.5, 1.2, 1.8),
  series_names = "my_series")

}
\seealso{
[iAR], [CiAR]
}
