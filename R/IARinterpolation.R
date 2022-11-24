#' Interpolation from IAR model
#'
#' Interpolation of missing values from models fitted by \code{\link{IARkalman}}
#'
#' @param x A given phi coefficient of the IAR model.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#' @param zero.mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the IAR process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{IARkalman}}
#'
#'
#' @examples
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' phi=IARkalman(y=y,st=st)$phi
#' print(phi)
#' napos=10
#' y0=y
#' y[napos]=NA
#' xest=phi
#' yest=IARinterpolation(xest,y=y,st=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted,col="red",pch=20)
IARinterpolation<-function(x,y,st,delta=0,yini=0,zero.mean=TRUE,standardized=TRUE)
{
  aux<-1e10
  value<-1e10
  br<-0
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  if (yini==0)
    yini = rnorm(1)
  out = optim(c(yini), IARphikalman, lower = -Inf, upper = Inf,
              x=x,y = y, st = st, yerr=delta,zeroMean = zero.mean,standardized=standardized, method="BFGS")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
