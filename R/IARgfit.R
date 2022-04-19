#' Fitted Values of IAR-Gamma model
#'
#' Fit an IAR-Gamma model to an irregularly observed time series.
#'
#' @param phi Estimated phi parameter by the iAR-Gamma model.
#' @param mu Estimated mu parameter by the iAR-Gamma model.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#'
#' @return Fitted values of the iAR-Gamma model
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{IARgsample}}, \code{\link{IARgamma}}
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
#' model<-IARgamma(y$y, st=st)
#' phi=model$phi
#' muest=model$mu
#' sigmaest=model$sigma
#' fit=IARgfit(phi=phi,mu=muest,y=y$y,st=st)
IARgfit<-function(phi,mu,y,st)
{
  delta=diff(st)
  y1=y[-c(length(y))]
  xd=phi**delta
  yhat = mu + xd * y1
  fit=c(mu,yhat)
  return(fitted=fit)
}
