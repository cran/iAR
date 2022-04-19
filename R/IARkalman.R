#' Maximum Likelihood Estimation of the IAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the IAR model parameter phi. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{arima}},\code{\link{IARphikalman}}
#'
#' @examples
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' phi=IARkalman(y=y,st=st)$phi
#' print(phi)
IARkalman<-function (y, st, delta=0,zero.mean='FALSE',standarized='TRUE')
{
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  out = optimize(IARphikalman, interval = c(0, 1), y = y, st = st, yerr=delta,zeroMean = zero.mean,standarized=standarized)
  phi = out$minimum
  ll = out$objective
  return(list(phi = phi, kalman = ll))
}
