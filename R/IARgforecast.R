#' Forecast from IAR-Gamma model
#'
#' Forecast from models fitted by \code{\link{IARgamma}}
#'
#' @param phi Estimated phi parameter by the iAR-Gamma model.
#' @param mu Estimated mu parameter by the iAR-Gamma model.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param tAhead The time ahead for forecast is required.
#'
#' @return Forecasted value from the iAR-Gamma model
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{IARgsample}}, \code{\link{IARgamma}}, \code{\link{IARgfit}}
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
#' y<-y$y
#' n=length(y)
#' p=trunc(n*0.99)
#' ytr=y[1:p]
#' yte=y[(p+1):n]
#' str=st[1:p]
#' ste=st[(p+1):n]
#' tahead=ste-str[p]
#' model<-IARgamma(ytr, st=str)
#' phi=model$phi
#' muest=model$mu
#' sigmaest=model$sigma
#' fit=IARgforecast(phi=phi,mu=muest,y=ytr,st=str,tAhead=tahead)
IARgforecast<-function(phi,mu,y,st,tAhead)
{
  y1=y[c(length(y))]
  xd=phi**(tAhead)
  yhat = mu + xd * y1
  fit=c(yhat)
  return(fitted=fit)
}
