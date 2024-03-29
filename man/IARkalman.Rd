% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IARkalman.R
\name{IARkalman}
\alias{IARkalman}
\title{Maximum Likelihood Estimation of the IAR Model via Kalman Recursions}
\usage{
IARkalman(y, st, delta = 0, zero.mean = TRUE, standardized = TRUE)
}
\arguments{
\item{y}{Array with the time series observations.}

\item{st}{Array with the irregular observational times.}

\item{delta}{Array with the measurements error standard deviations.}

\item{zero.mean}{logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.}

\item{standardized}{logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.}
}
\value{
A list with the following components:
\itemize{
\item{phi}{ MLE of the phi parameter of the IAR model.}
\item{ll}{ Value of the negative log likelihood evaluated in phi.}
}
}
\description{
Maximum Likelihood Estimation of the IAR model parameter phi. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
}
\examples{
set.seed(6714)
st<-gentime(n=100)
y<-IARsample(phi=0.99,st=st,n=100)
y<-y$series
phi=IARkalman(y=y,st=st)$phi
print(phi)
}
\references{
\insertRef{Eyheramendy_2018}{iAR}
}
\seealso{
\code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{arima}},\code{\link{IARphikalman}}
}
