#' Maximum Likelihood Estimation of the CIAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the CIAR model parameters phi.R and phi.I. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#' @param niter Number of iterations in which the function nlminb will be repeated.
#' @param seed a single value, interpreted as the seed of the random process.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phiR}{ MLE of the Real part of the coefficient of CIAR model (phi.R).}
#' \item{phiI}{ MLE of the Imaginary part of the coefficient of the CIAR model (phi.I).}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi.R and phi.I.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIAR.sample}}, \code{\link{CIAR.phi.kalman}}
#'
#'
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIAR.sample(n=n,phi.R=0.9,phi.I=0,sT=st,c=1)
#' y=x$y
#' y1=y/sd(y)
#' ciar=CIAR.kalman(y=y1,t=st)
#' ciar
#' Mod(complex(real=ciar$phiR,imaginary=ciar$phiI))
CIAR.kalman<-function(y,t,delta=0,zero.mean='TRUE',standarized='TRUE',c=1,niter=10,seed=1234)
{
  set.seed(seed)
  aux<-1e10
  value<-1e10
  br<-0
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  for(i in 1:niter)
  {
    phi.R=2*runif(1)-1
    phi.I=2*runif(1)-1
    #print(complex(1,real=phi.R,imaginary=phi.I))
    if(Mod(complex(1,real=phi.R,imaginary=phi.I))<1)
    {
      optim<-nlminb(start=c(phi.R,phi.I),objective=CIAR.phi.kalman,y=y,t=t,yerr=delta,zero.mean=zero.mean,standarized=standarized,c=c,lower=c(-1,-1),upper=c(1,1))
      value<-optim$objective
    }
    if(aux>value)
    {
      par<-optim$par
      aux<-value
      br<-br+1
    }
    if(aux<=value & br>1 & i>trunc(niter/2))
      break;
  }
  if(aux==1e10)
    par<-c(0,0)
  return(list(phiR=par[1],phiI=par[2],ll=aux))
}
