#' Minus Log Likelihood IAR-Gamma Model
#'
#' This function return the negative log likelihood of the IAR-Gamma given specific values of phi, mu and sigma.
#'
#' @param x An array with the parameters of the IAR-Gamma model. The first element of the array corresponding to the phi parameter, the second to the level parameter mu, and the last one to the scale parameter sigma.
#' @param y Array with the time series observations.
#' @param sT Array with the irregular observational times.
#'
#' @return Value of the negative log likelihood evaluated in phi, mu and sigma.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARg.sample}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARg.sample(n,phi=0.9,st,sigma2=1,mu=1)
#' IAR.phi.gamma(x=c(0.9,1,1),y=y$y,sT=st)
IAR.phi.gamma<-function (x, y, sT)
{
  mu=x[2]
  sigma=x[3]
  x=x[1]
  n = length(y)
  d <- diff(sT)
  xd=x**d
  yhat = mu+xd * y[-n]  #Mean of conditional distribution
  gL=sigma*(1-xd**(2))  #Variance of conditional distribution
  beta=gL/yhat #Beta parameter of gamma distribution
  alpha=yhat**2/gL #Alpha parameter of gamma distribution
  out=sum((-alpha)*log(beta) - lgamma(alpha) - y[-1]/beta + (alpha-1)*log(y[-1])) - y[1]  #Full Log Likelihood
  out=-out #-Log Likelihood (We want to minimize it)
  return(out)
}
