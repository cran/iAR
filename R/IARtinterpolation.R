#' Interpolation from IAR-T model
#'
#' Interpolation of missing values from models fitted by \code{\link{IARt}}
#'
#' @param x A given array with the parameters of the IAR-T model. The first element of the array corresponding to the phi parameter and the second element to the scale parameter sigma
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param nu degrees of freedom
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the IAR-T process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARtsample}}, \code{\link{IARt}}
#'
#'
#' @examples
#' set.seed(6714)
#' n<-100
#' st<-gentime(n)
#' y<-IARtsample(n,0.9,st,sigma2=1,nu=3)
#' model<-IARt(y$y, st=st)
#' napos=10
#' y0=y$y
#' y=y$y
#' y[napos]=NA
#' xest=c(model$phi,model$sigma)
#' yest=IARtinterpolation(x=xest,y=y,st=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted,col="red",pch=20)
IARtinterpolation<-function(x,y,st,nu=3,yini=0)
{
  aux<-1e10
  value<-1e10
  br<-0
  if (yini==0)
    yini = rnorm(1)
  out = optim(yini, IARphit, lower = -Inf, upper = Inf,
              x=x,y=y,st=st,nu=nu, method="L-BFGS-B")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
