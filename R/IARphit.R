#' Minus Log Likelihood IAR-T Model
#'
#' This function return the negative log likelihood of the IAR-T given specific values of phi and sigma.
#'
#' @param x An array with the parameters of the IAR-T model. The first element of the array corresponding to the phi parameter and the second element to the scale parameter sigma
#' @param y Array with the time series observations
#' @param sT Array with the irregular observational times
#' @param nu degrees of freedom
#'
#' @return Value of the negative log likelihood evaluated in phi,sigma and nu.
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
#' st<-gentime(n) #Unequally spaced times
#' y<-IARt.sample(n,0.9,st,sigma2=1,nu=3)
#' IAR.phi.t(x=c(0.9,1),y=y$y,sT=st)
IAR.phi.t<-function (x, y, sT, nu=3) #Minus Log Full Likelihood Function
{
  sigma=x[2]
  x=x[1]
  n = length(y)
  d <- diff(sT)
  xd=x**d
  yhat = xd * y[-n]  #Mean of conditional distribution
  gL=sigma*(1-xd**(2))*((nu-2)/nu)  #Variance of conditional distribution
  cte = (n-1)*log((gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi))))
  stand=((y[-1]-yhat)/sqrt(gL))**2
  s1=sum(0.5*log(gL))
  s2=sum(log(1 + (1/nu)*stand))
  out= cte - s1 - ((nu+1)/2)*s2 -0.5*(log(2*pi) + y[1]**2)
  out=-out #-Log Likelihood (We want to minimize it)
  return(out)
}
