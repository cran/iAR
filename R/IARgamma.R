#' Maximum Likelihood Estimation of the IAR-Gamma model
#'
#' Maximum Likelihood Estimation of the IAR-Gamma model.
#'
#' @param y Array with the time series observations
#' @param st Array with the irregular observational times
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR-Gamma model.}
#' \item{mu}{ MLE of the mu parameter of the IAR-Gamma model.}
#' \item{sigma}{ MLE of the sigma parameter of the IAR-Gamma model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi, mu and sigma.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARgsample}}, \code{\link{IARphigamma}}
#'
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
#' model<-IARgamma(y$y, st=st)
#' phi=model$phi
#' muest=model$mu
#' sigmaest=model$sigma
IARgamma<-function(y, st)
{
  aux<-1e10
  value<-1e10
  br<-0
  for(i in 1:20)
  {
    phi=runif(1)
    mu=mean(y)*runif(1)
    sigma=var(y)*runif(1)
    optim<-nlminb(start=c(phi,mu,sigma),objective=IARphigamma,y=y,st=st,lower=c(0,0.0001,0.0001),upper=c(0.9999,mean(y),var(y)))
    value<-optim$objective
    if(aux>value)
    {
      par<-optim$par
      aux<-value
      br<-br+1
    }
    if(aux<=value & br>5 & i>10)
      break;
  }
  if(aux==1e10)
    par<-c(0,0,0)
  return(list(phi=par[1],mu=par[2],sigma=par[3],ll=aux))
}
