#' Simulate from an IAR-T Model
#'
#' Simulates an IAR-T Time Series Model.
#'
#' @param n Length of the output time series. A strictly positive integer.
#' @param phi A coefficient of IAR-T model. A value between 0 and 1.
#' @param st Array with observational times.
#' @param sigma2 Scale parameter of the IAR-T process. A positive value.
#' @param nu degrees of freedom.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{y}{ Array with simulated IAR-t process.}
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
#' y<-IARtsample(n,0.9,st,sigma2=1,nu=3)
#' plot(st,y$y,type='l')
#' hist(y$y,breaks=20)
IARtsample<-function(n,phi,st,sigma2=1,nu=3)
{
  delta<-diff(st) #Times differences
  y <- vector("numeric", n)
  y[1] <- rnorm(1) #initialization
  for (i in 2:n)
  {
    phid=phi**(delta[i-1]) #Expected Value Conditional Distribution
    yhat = phid * y[i-1]  #Mean of conditional distribution
    gL=sigma2*(1-phid**(2)) #Variance Value Conditional Distribution
    y[i] <- rt(1, df=nu)*sqrt(gL * (nu-2)/nu) + yhat #Conditional-t IAR
  }
  return(list(y=y,st=st))
}
