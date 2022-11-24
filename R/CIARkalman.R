#' Maximum Likelihood Estimation of the CIAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the CIAR model parameters phiR and phiI. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param zero.mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#' @param niter Number of iterations in which the function nlminb will be repeated.
#' @param seed a single value, interpreted as the seed of the random process.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phiR}{ MLE of the Real part of the coefficient of CIAR model (phiR).}
#' \item{phiI}{ MLE of the Imaginary part of the coefficient of the CIAR model (phiI).}
#' \item{ll}{ Value of the negative log likelihood evaluated in phiR and phiI.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIARsample}}, \code{\link{CIARphikalman}}
#'
#'
#' @examples
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
#' y=x$y
#' y1=y/sd(y)
#' ciar=CIARkalman(y=y1,t=st)
#' ciar
#' Mod(complex(real=ciar$phiR,imaginary=ciar$phiI))
CIARkalman<-function(y,t,delta=0,zero.mean=TRUE,standardized=TRUE,c=1,niter=10,seed=1234)
{
  set.seed(seed)
  aux<-1e10
  value<-1e10
  br<-0
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  for(i in 1:niter)
  {
    phiR=2*runif(1)-1
    phiI=2*runif(1)-1
    if(Mod(complex(1,real=phiR,imaginary=phiI))<1)
    {
      optim <- nlminb(start=c(phiR,phiI), objective=CIARphikalman, y=y, t=t,
                      yerr=delta, zeroMean=zero.mean, standardized=standardized,
                      c=c,yest=0, lower=c(-1,-1), upper=c(1,1))

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
