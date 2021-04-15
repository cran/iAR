#' Simulate from an IAR-Gamma Model
#'
#' Simulates an IAR-Gamma Time Series Model.
#'
#' @param n Length of the output time series. A strictly positive integer.
#' @param phi A coefficient of IAR-Gamma model. A value between 0 and 1.
#' @param st Array with observational times.
#' @param sigma2 Scale parameter of the IAR-Gamma process. A positive value.
#' @param mu Level parameter of the IAR-Gamma process. A positive value.
#'
#' @return  A list with the following components:
#' \itemize{
#' \item{y}{ Array with simulated IAR-Gamma process.}
#' \item{st}{ Array with observation times.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARg.sample(n,phi=0.9,st,sigma2=1,mu=1)
#' plot(st,y$y,type='l')
#' hist(y$y,breaks=20)
#'
IARg.sample<-function(n,phi,st,sigma2=1,mu=1)
{
  delta<-diff(st) #Times differences
  y <- vector("numeric", n)
  y[1] <- rgamma(1,shape=1,scale=1) #initialization
  shape<-rep(0,n)
  scale<-rep(0,n)
  yhat<-rep(0,n)
  for (i in 2:n)
  {
    phid=phi**(delta[i-1]) #Expected Value Conditional Distribution
    yhat[i] = mu+phid * y[i-1]  #Mean of conditional distribution
    gL=sigma2*(1-phid**(2)) #Variance Value Conditional Distribution
    shape[i]=yhat[i]**2/gL
    scale[i]=(gL/yhat[i])
    y[i] <- rgamma(1,shape=shape[i], scale=scale[i])#Conditional Gamma IAR
  }
  #ts.plot(yhat)
  return(list(y=y,st=st))
}
