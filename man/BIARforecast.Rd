% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{BIARforecast}
\alias{BIARforecast}
\title{Forecast from BIAR model}
\usage{
BIARforecast(phiR, phiI, y1, y2, t, tAhead)
}
\arguments{
\item{phiR}{Autocorrelation coefficient of BIAR model.}

\item{phiI}{Cross-correlation coefficient of BIAR model.}

\item{y1}{Array with the observations of the first time series of the BIAR process.}

\item{y2}{Array with the observations of the second time series of the BIAR process.}

\item{t}{Array with the observational times.}

\item{tAhead}{The time ahead for which the forecast is required.}
}
\value{
A list with the following components:
\itemize{
\item{fitted}{ Fitted values by the BIAR model.}
\item{forecast}{ Point forecast in the time ahead required.}
\item{Lambda}{ Lambda value estimated by the BIAR model at the last time point.}
\item{Sighat}{ Covariance matrix estimated by the BIAR model at the last time point.}
}
}
\description{
Forecast from models fitted by \code{\link{BIARkalman}}
}
\examples{
#Simulated Data
n=100
set.seed(6714)
st<-gentime(n)
x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
biar=iAR::BIARkalman(y1=x$y[1,],y2=x$y[2,],t=st)
forBIAR<-BIARforecast(phiR=biar$phiR,phiI=biar$phiI,y1=x$y[1,],y2=x$y[2,],t=st,tAhead=c(1.3))
}
\references{
\insertRef{Elorrieta_2021}{iAR}
}
\seealso{
\code{\link{BIARsample}}, \code{\link{BIARkalman}}, \code{\link{BIARfit}}
}
