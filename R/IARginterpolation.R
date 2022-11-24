#' Interpolation from IAR-Gamma model
#'
#' Interpolation of missing values from models fitted by \code{\link{IARgamma}}
#'
#' @param x A given array with the parameters of the IAR-Gamma model. The first element of the array corresponding to the phi parameter, the second to the level parameter mu, and the last one to the scale parameter sigma.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the IAR-Gamma process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARgsample}}, \code{\link{IARgamma}}
#'
#'
#' @examples
#' set.seed(6714)
#' n<-100
#' st<-gentime(n)
#' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
#' model<-IARgamma(y$y, st=st)
#' y<-y$y
#' napos=10
#' y0=y
#' y[napos]=NA
#' xest=c(model$phi,model$mu,model$sigma)
#' yest=IARginterpolation(x=xest,y=y,st=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted,col="red",pch=20)
IARginterpolation<-function(x,y,st,yini=1)
{
  aux<-1e10
  value<-1e10
  br<-0
  if (yini==1)
    yini = rgamma(1,1,1)
  out = optim(yini, IARphigamma, lower = 0.0001, upper = Inf,
              x_input=x,y = y, st = st, method="L-BFGS-B")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
