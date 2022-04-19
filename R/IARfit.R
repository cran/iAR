#' Fitted Values of IAR model
#'
#' Fit an IAR model to an irregularly observed time series.
#'
#' @param phi Estimated phi parameter by the iAR model.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series
#' @param include.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#'
#' @return Fitted values of the iAR model
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{IARloglik}}, \code{\link{IARkalman}}
#' @examples
#'
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' phi=IARloglik(y=y,st=st)$phi
#' fit=IARfit(phi=phi,y=y,st=st)
IARfit<-function(phi,y,st,standarized=T,include.mean=F)
{
  sigma = 1
  mu = 0
  if (standarized == F)
    sigma = var(y)
  if (include.mean == T)
    mu = mean(y)
  if (include.mean == T)
    y = y-mu
  if (standarized == F)
    y=y/sqrt(sigma)
  delta=diff(st)
  y1=y[-c(length(y))]
  fit=c(0,(phi**delta)*y1)
  fit=(fit*sqrt(sigma)+mu)
  return(fitted=fit)
}

