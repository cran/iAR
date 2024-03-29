% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{IARphikalman}
\alias{IARphikalman}
\title{Minus Log Likelihood of the IAR Model estimated via Kalman Recursions}
\usage{
IARphikalman(yest, x, y, yerr, st, zeroMean = TRUE, standardized = TRUE)
}
\arguments{
\item{yest}{The estimate of a missing value in the time series. This function recognizes a missing value with a NA. If the time series does not have a missing value, this value does not affect the computation of the likelihood.}

\item{x}{A given phi coefficient of the IAR model.}

\item{y}{Array with the time series observations.}

\item{yerr}{Array with the measurements error standard deviations.}

\item{st}{Array with the irregular observational times.}

\item{zeroMean}{logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.}

\item{standardized}{logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.}
}
\value{
Value of the negative log likelihood evaluated in phi.
}
\description{
This function return the negative log likelihood of the IAR process given a specific value of phi.
}
\examples{
set.seed(6714)
st<-gentime(n=100)
y<-IARsample(phi=0.99,st=st,n=100)
y<-y$series
yerr=rep(0,100)
IARphikalman(x=0.8,y=y,yerr=yerr,st=st,yest=0)
}
\references{
\insertRef{Eyheramendy_2018}{iAR}
}
\seealso{
\code{\link{gentime}}, \code{\link{IARsample}}
}
