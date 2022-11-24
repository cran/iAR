#' Interpolation from CIAR model
#'
#' Interpolation of missing values from models fitted by \code{\link{CIARkalman}}
#'
#' @param x An array with the parameters of the CIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of CIAR model.
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#' @param zero.mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#' @param seed a single value, interpreted as the seed of the random process.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the CIAR process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIARsample}}, \code{\link{CIARkalman}}
#'
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
#' napos=10
#' y0=y1
#' y1[napos]=NA
#' xest=c(ciar$phiR,ciar$phiI)
#' yest=CIARinterpolation(xest,y=y1,t=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted*sd(y),col="red",pch=20)
CIARinterpolation<-function(x,y,t,delta=0,yini=0,zero.mean=TRUE,standardized=TRUE,c=1,seed=1234)
{
  set.seed(seed)
  aux<-1e10
  value<-1e10
  br<-0
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  if (yini==0)
    yini = rnorm(1)
  out = optim(c(yini), CIARphikalman, lower = -Inf, upper = Inf,
              x=x, y=y, t=t,yerr=delta, zeroMean=zero.mean, standardized=standardized,
              c=c, method="BFGS")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
