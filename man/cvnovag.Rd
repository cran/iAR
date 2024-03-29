% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvnovag.R
\docType{data}
\name{cvnovag}
\alias{cvnovag}
\title{ZTF g-band Cataclysmic Variable/Nova}
\format{
A data frame with 67 observations on the following 3 variables:
\describe{
  \item{t}{heliocentric Julian Day - 2400000}
  \item{m}{magnitude}
  \item{merr}{measurement error standard deviations.}
}
}
\usage{
cvnovag
}
\description{
Time series of a cataclysmic variable/nova object observed in the g-band of the ZTF survey and processed by the ALeRCE broker.ZTF Object code: ZTF18aayzpbr
}
\examples{
data(cvnovag)
plot(cvnovag$t,cvnovag$m,type="l",ylab="",xlab="",col="green")
}
\references{
\insertRef{Forster_2021}{iAR}
}
\keyword{datasets}
