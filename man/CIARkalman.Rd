% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIARkalman.R
\name{CIARkalman}
\alias{CIARkalman}
\title{Maximum Likelihood Estimation of the CIAR Model via Kalman Recursions}
\usage{
CIARkalman(
  y,
  t,
  delta = 0,
  zero.mean = TRUE,
  standardized = TRUE,
  c = 1,
  niter = 10,
  seed = 1234
)
}
\arguments{
\item{y}{Array with the time series observations.}

\item{t}{Array with the irregular observational times.}

\item{delta}{Array with the measurements error standard deviations.}

\item{zero.mean}{logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.}

\item{standardized}{logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.}

\item{c}{Nuisance parameter corresponding to the variance of the imaginary part.}

\item{niter}{Number of iterations in which the function nlminb will be repeated.}

\item{seed}{a single value, interpreted as the seed of the random process.}
}
\value{
A list with the following components:
\itemize{
\item{phiR}{ MLE of the Real part of the coefficient of CIAR model (phiR).}
\item{phiI}{ MLE of the Imaginary part of the coefficient of the CIAR model (phiI).}
\item{ll}{ Value of the negative log likelihood evaluated in phiR and phiI.}
}
}
\description{
Maximum Likelihood Estimation of the CIAR model parameters phiR and phiI. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
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
Mod(complex(real=ciar$phiR,imaginary=ciar$phiI))
}
\references{
\insertRef{Elorrieta_2019}{iAR}
}
\seealso{
\code{\link{gentime}}, \code{\link{CIARsample}}, \code{\link{CIARphikalman}}
}
