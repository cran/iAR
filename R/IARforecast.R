#' Forecast from IAR model
#'
#' Forecast from models fitted by \code{\link{IARloglik}}
#'
#' @param phi Estimated phi parameter by the iAR model.
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series
#' @param zero.mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param tAhead The time ahead for forecast is required.
#'
#' @return Forecasted value from the iAR model
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{IARloglik}}, \code{\link{IARkalman}}, \code{\link{IARfit}}
#' @examples
#'
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' n=length(y)
#' p=trunc(n*0.99)
#' ytr=y[1:p]
#' yte=y[(p+1):n]
#' str=st[1:p]
#' ste=st[(p+1):n]
#' tahead=ste-str[p]
#' phi=IARloglik(y=ytr,st=str)$phi
#' forIAR=IARforecast(phi=phi,y=ytr,st=str,tAhead=tahead)
IARforecast<-function(phi,y,st,standardized=TRUE,zero.mean=TRUE,tAhead)
{
  sigma = 1
  mu = 0
  if (standardized == FALSE)
    sigma = var(y)
  if (zero.mean == FALSE)
    mu = mean(y)
  if (zero.mean == FALSE)
    y = y-mu
  if (standardized == FALSE)
    y=y/sqrt(sigma)
  y1=y[length(y)]
  fit=(phi**(tAhead)*y1)
  fit=(fit*sqrt(sigma)+mu)
  return(fitted=fit)
}

