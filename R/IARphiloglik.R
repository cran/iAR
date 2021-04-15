#' Minus Log Likelihood of the IAR Model
#'
#' This function return the negative log likelihood of the IAR Model for a specific value of phi.
#'
#' @param x A given phi coefficient of the IAR model.
#' @param y Array with the time series observations.
#' @param sT Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param include.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y was standarized; if false, y contains the raw data
#'
#' @return Value of the negative log likelihood evaluated in phi.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IAR.sample}}
#'
#' @examples
#'
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IAR.sample(phi=0.99,n=100,st)
#' y<-y$series
#' IAR.phi.loglik(x=0.8,y=y,sT=st)
IAR.phi.loglik=function(x,y,sT,delta=0,include.mean='FALSE',standarized='TRUE')
{
  sigma=1
  mu=0
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  if(standarized=='FALSE')
    sigma=var(y)
  if(include.mean=='TRUE')
    mu=mean(y)
  n=length(y)
  delta<-delta[-1]
  d<-diff(sT)
  phi=x**d
  yhat=mu+phi*(y[-n]-mu)
  cte=(n/2)*log(2*pi)
  s1=cte+0.5*sum(log(sigma*(1-phi**2)+delta**2)+(y[-1]-yhat)**2/(sigma*(1-phi**2)+delta**2))
  return(s1)
}
