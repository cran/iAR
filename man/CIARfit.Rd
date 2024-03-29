% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CIARfit}
\alias{CIARfit}
\title{Fitted Values of CIAR model}
\usage{
CIARfit(phiValues, y, t, standardized = TRUE, c = 1)
}
\arguments{
\item{phiValues}{An array with the parameters of the CIAR model. The elements of the array are, in order, the real and the imaginary part of the phi parameter of the CIAR model.}

\item{y}{Array with the time series observations.}

\item{t}{Array with the irregular observational times.}

\item{standardized}{logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series}

\item{c}{Nuisance parameter corresponding to the variance of the imaginary part.}
}
\value{
A list with the following components:
\itemize{
\item{yhat}{ Fitted values of the observable part of CIAR model.}
\item{xhat}{ Fitted values of both observable part and imaginary part of CIAR model.}
\item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
\item{Theta}{ Theta array estimated by the CIAR model at the last time point.}
\item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
\item{Qt}{ Covariance matrix of the state equation estimated by the CIAR model at the last time point.}
}
}
\description{
Fit a CIAR model to an irregularly observed time series.
}
\examples{
n=100
set.seed(6714)
st<-gentime(n)
x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
y=x$y
y1=y/sd(y)
ciar=CIARkalman(y=y1,t=st)
ciar
yhat=CIARfit(phiValues=c(ciar$phiR,ciar$phiI),y=y1,t=st)
}
\references{
\insertRef{Elorrieta_2019}{iAR}
}
\seealso{
\code{\link{gentime}}, \code{\link{CIARsample}}, \code{\link{CIARphikalman}},\code{\link{CIARkalman}}
}
