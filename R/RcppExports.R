# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Full Minus Log Likelihood of the BIAR Model
#'
#' This function return the full negative log likelihood of the BIAR process given specific values of phiR and phiI
#'
#'
#' @param yest An array with the estimate of a missing value in one or both time series of the bivariate process. This function recognizes a missing value with a NA. If the bivariate time series does not have a missing value, this value does not affect the computation of the likelihood.
#' @param phiValues An array with the parameters of the BIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of BIAR model.
#' @param y1 Array with the observations of the first time series of the BIAR process.
#' @param y2 Array with the observations of the second time series of the BIAR process.
#' @param t Array with the irregular observational times.
#' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
#' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
#' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#'
#' @return Value of the full negative log likelihood evaluated in phiR and phiI.
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{BIARsample}}
#'
#'
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
#' y=x$y
#' y1=y[1,]
#' y2=y[2,]
#' yerr1=rep(0,n)
#' yerr2=rep(0,n)
#' BIARLL(yest=0,phiValues=c(0.8,0.2),y1=y1,y2=y2,t=st,yerr1=yerr1,yerr2=yerr2)
BIARLL <- function(yest, phiValues, y1, y2, t, yerr1, yerr2, zeroMean = "TRUE") {
    .Call(`_iAR_BIARLL`, yest, phiValues, y1, y2, t, yerr1, yerr2, zeroMean)
}

#' Fitted Values of BIAR model
#'
#' Fit a BIAR model to a bivariate irregularly observed time series.
#'
#' @param phiValues An array with the parameters of the BIAR model. The elements of the array are, in order, the autocorrelation and the cross correlation parameter of the BIAR model.
#' @param y1 Array with the observations of the first time series of the BIAR process.
#' @param y2 Array with the observations of the second time series of the BIAR process.
#' @param t Array with the irregular observational times.
#' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
#' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
#' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{rho}{ Estimated value of the contemporary correlation coefficient.}
#' \item{innov.var}{ Estimated value of the innovation variance.}
#' \item{fitted}{ Fitted values of the BIAR model.}
#' \item{fitted.state}{ Fitted state values of the BIAR model.}
#' \item{Lambda}{ Lambda value estimated by the BIAR model at the last time point.}
#' \item{Theta}{ Theta array estimated by the BIAR model at the last time point.}
#' \item{Sighat}{ Covariance matrix estimated by the BIAR model at the last time point.}
#' \item{Qt}{ Covariance matrix of the state equation estimated by the BIAR model at the last time point.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#' @seealso
#' \code{\link{gentime}}, \code{\link{BIARsample}}, \code{\link{BIARphikalman}}, \code{\link{BIARkalman}}
#'
#' @examples
#' \donttest{
#' n=80
#' set.seed(6714)
#' st<-gentime(n)
#' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st,rho=0.9)
#' y=x$y
#' y1=y/apply(y,1,sd)
#' yerr1=rep(0,n)
#' yerr2=rep(0,n)
#' biar=BIARkalman(y1=y1[1,],y2=y1[2,],t=st,delta1 = yerr1,delta2=yerr2)
#' biar
#' predbiar=BIARfit(phiValues=c(biar$phiR,biar$phiI),y1=y1[1,],y2=y1[2,],t=st,yerr1
#'  = rep(0,length(y[1,])),yerr2=rep(0,length(y[1,])))
#' rho=predbiar$rho
#' print(rho)
#' yhat=predbiar$fitted
#' }
BIARfit <- function(phiValues, y1, y2, t, yerr1, yerr2, zeroMean = "TRUE") {
    .Call(`_iAR_BIARfit`, phiValues, y1, y2, t, yerr1, yerr2, zeroMean)
}

#' Minus Log Likelihood of the BIAR Model
#'
#' This function return the negative log likelihood of the BIAR process given specific values of phiR and phiI
#'
#'
#' @param phiValues An array with the parameters of the BIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of BIAR model.
#' @param y1 Array with the observations of the first time series of the BIAR process.
#' @param y2 Array with the observations of the second time series of the BIAR process.
#' @param t Array with the irregular observational times.
#' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
#' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
#' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#'
#' @return Value of the negative log likelihood evaluated in phiR and phiI.
#' @export
#'
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#' @seealso
#' \code{\link{gentime}}, \code{\link{BIARsample}}
#'
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
#' y=x$y
#' y1=y[1,]
#' y2=y[2,]
#' yerr1=rep(0,n)
#' yerr2=rep(0,n)
#' BIARphikalman(phiValues=c(0.8,0.2),y1=y1,y2=y2,t=st,yerr1=yerr1,yerr2=yerr2)
BIARphikalman <- function(phiValues, y1, y2, t, yerr1, yerr2, zeroMean = "TRUE") {
    .Call(`_iAR_BIARphikalman`, phiValues, y1, y2, t, yerr1, yerr2, zeroMean)
}

#' Fitted Values of CIAR model
#'
#' Fit a CIAR model to an irregularly observed time series.
#'
#' @param phiValues An array with the parameters of the CIAR model. The elements of the array are, in order, the real and the imaginary part of the phi parameter of the CIAR model.
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{yhat}{ Fitted values of the observable part of CIAR model.}
#' \item{xhat}{ Fitted values of both observable part and imaginary part of CIAR model.}
#' \item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
#' \item{Theta}{ Theta array estimated by the CIAR model at the last time point.}
#' \item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
#' \item{Qt}{ Covariance matrix of the state equation estimated by the CIAR model at the last time point.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIARsample}}, \code{\link{CIARphikalman}},\code{\link{CIARkalman}}
#'
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
#' y=x$y
#' y1=y/sd(y)
#' ciar=CIARkalman(y=y1,t=st)
#' ciar
#' yhat=CIARfit(phiValues=c(ciar$phiR,ciar$phiI),y=y1,t=st)
CIARfit <- function(phiValues, y, t, standarized = "TRUE", c = 1) {
    .Call(`_iAR_CIARfit`, phiValues, y, t, standarized, c)
}

#' Forecast from CIAR model
#'
#' Forecast from models fitted by \code{\link{CIARkalman}}
#'
#' @param phiR Real part of the phi coefficient of CIAR model.
#' @param phiI Imaginary part of the phi coefficient of CIAR model.
#' @param y1 Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param nAhead The number of steps ahead for forecast is required.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Fitted values by the CIAR model.}
#' \item{forecast}{ Point Forecasts in the n.ahead times.}
#' \item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
#' \item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{CIARsample}}, \code{\link{CIARkalman}}, \code{\link{CIARfit}}
#'
#'
#' @examples
#' #Simulated Data
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
#' y=x$y
#' y1=y/sd(y)
#' n=length(y1)
#' p=trunc(n*0.99)
#' ytr=y1[1:p]
#' yte=y1[(p+1):n]
#' str=st[1:p]
#' ste=st[(p+1):n]
#' n.ahead=ste-str[p]
#'
#' final<-matrix(0,length(n.ahead),4)
#' ciar=CIARkalman(y=ytr,t=str)
#' forCIAR<-CIARforecast(ciar$phiR,ciar$phiI,ytr,str,nAhead=n.ahead)
CIARforecast <- function(phiR, phiI, y1, st, nAhead) {
    .Call(`_iAR_CIARforecast`, phiR, phiI, y1, st, nAhead)
}

#' Minus Log Likelihood of the CIAR Model
#'
#' This function return the negative log likelihood of the CIAR process given specific values of phiR and phiI
#'
#' @param x An array with the parameters of the CIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of CIAR model.
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param yerr Array with the measurements error standard deviations.
#' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @return Value of the negative log likelihood evaluated in phiR and phiI.
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIARsample}}
#'
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
#' y=x$y
#' yerr=rep(0,n)
#' CIARphikalman(x=c(0.8,0),y=y,t=st,yerr=yerr)
CIARphikalman <- function(x, y, t, yerr, zeroMean = "TRUE", standarized = "TRUE", c = 1.0) {
    .Call(`_iAR_CIARphikalman`, x, y, t, yerr, zeroMean, standarized, c)
}

#' Simulate from a CIAR Model
#'
#' Simulates a CIAR Time Series Model
#'
#' @param n Length of the output time series. A strictly positive integer.
#' @param st Array with observational times.
#' @param phiR Real part of the coefficient of CIAR model. A value between -1 and 1.
#' @param phiI Imaginary part of the coefficient of CIAR model. A value between -1 and 1.
#' @param rho Correlation between the real and the imaginary part of the process. A value between -1 and 1.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @details The chosen phiR and phiI values must satisfy the condition $|phiR + i phiI| < 1$.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{y}{Array with the simulated real part of the CIAR process.}
#' \item{t}{ Array with observation times.}
#' \item{Sigma}{ Covariance matrix of the process.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
#' plot(st,x$y,type='l')
#' x=CIARsample(n=n,phiR=-0.9,phiI=0,st=st,c=1)
#' plot(st,x$y,type='l')
CIARsample <- function(n, phiR, phiI, st, rho = 0L, c = 1L) {
    .Call(`_iAR_CIARsample`, n, phiR, phiI, st, rho, c)
}

#' Simulate from an IAR-Gamma Model
#'
#' Simulates an IAR-Gamma Time Series Model.
#'
#' @param phi A coefficient of IAR-Gamma model. A value between 0 and 1.
#' @param st Array with observational times.
#' @param n Length of the output time series. A strictly positive integer.
#' @param sigma2 Scale parameter of the IAR-Gamma process. A positive value.
#' @param mu Level parameter of the IAR-Gamma process. A positive value.
#'
#' @return  A list with the following components:
#' \itemize{
#' \item{y}{ Array with simulated IAR-Gamma process.}
#' \item{st}{ Array with observation times.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
#' plot(st,y$y,type='l')
#' hist(y$y,breaks=20)
IARgsample <- function(phi, st, n = 100L, sigma2 = 1L, mu = 1L) {
    .Call(`_iAR_IARgsample`, phi, st, n, sigma2, mu)
}

#' Minus Log Likelihood IAR-Gamma Model
#'
#' This function return the negative log likelihood of the IAR-Gamma given specific values of phi, mu and sigma.
#'
#' @param x_input An array with the parameters of the IAR-Gamma model. The first element of the array corresponding to the phi parameter, the second to the level parameter mu, and the last one to the scale parameter sigma.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#'
#' @return Value of the negative log likelihood evaluated in phi, mu and sigma.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARgsample}}
#'
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
#' IARphigamma(x_input=c(0.9,1,1),y=y$y,st=st)
IARphigamma <- function(x_input, y, st) {
    .Call(`_iAR_IARphigamma`, x_input, y, st)
}

#' Minus Log Likelihood of the IAR Model estimated via Kalman Recursions
#'
#' This function return the negative log likelihood of the IAR process given a specific value of phi.
#'
#' @param x A given phi coefficient of the IAR model.
#' @param y Array with the time series observations.
#' @param yerr Array with the measurements error standard deviations.
#' @param st Array with the irregular observational times.
#' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
#'
#' @return Value of the negative log likelihood evaluated in phi.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARsample}}
#'
#' @examples
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' yerr=rep(0,100)
#' IARphikalman(x=0.8,y=y,yerr=yerr,st=st)
IARphikalman <- function(x, y, yerr, st, zeroMean = "TRUE", standarized = "TRUE") {
    .Call(`_iAR_IARphikalman`, x, y, yerr, st, zeroMean, standarized)
}

#' Minus Log Likelihood of the IAR Model
#'
#' This function return the negative log likelihood of the IAR Model for a specific value of phi.
#'
#' @param x A given phi coefficient of the IAR model.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param delta_input Array with the measurements error standard deviations.
#' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y was standarized; if false, y contains the raw data
#'
#' @return Value of the negative log likelihood evaluated in phi.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARsample}}
#'
#' @examples
#'
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' IARphiloglik(x=0.8,y=y,st=st,delta_input=c(0))
IARphiloglik <- function(x, y, st, delta_input, zeroMean = "TRUE", standarized = "TRUE") {
    .Call(`_iAR_IARphiloglik`, x, y, st, delta_input, zeroMean, standarized)
}

#' Minus Log Likelihood IAR-T Model
#'
#' This function return the negative log likelihood of the IAR-T given specific values of phi and sigma.
#'
#' @param x An array with the parameters of the IAR-T model. The first element of the array corresponding to the phi parameter and the second element to the scale parameter sigma
#' @param y Array with the time series observations
#' @param st Array with the irregular observational times
#' @param nu degrees of freedom
#'
#' @return Value of the negative log likelihood evaluated in phi,sigma and nu.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARgsample}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n) #Unequally spaced times
#' y<-IARtsample(n,0.9,st,sigma2=1,nu=3)
#' IARphit(x=c(0.9,1),y=y$y,st=st)
IARphit <- function(x, y, st, nu = 3) {
    .Call(`_iAR_IARphit`, x, y, st, nu)
}

#' Simulate from an IAR Model
#'
#' Simulates an IAR Time Series Model.
#'
#' @param phi A coefficient of IAR model. A value between 0 and 1
#' @param st Array with observational times.
#' @param n Length of the output time series. A strictly positive integer.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{times}{ Array with observation times.}
#' \item{series}{ Array with simulated IAR data.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @examples
#'
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st, n=100)
#' y<-y$series
#' plot(st,y,type='l')
IARsample <- function(phi, st, n = 100L) {
    .Call(`_iAR_IARsample`, phi, st, n)
}

