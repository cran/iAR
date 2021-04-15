#' Maximum Likelihood Estimation of the IAR-T model
#'
#' Maximum Likelihood Estimation of the IAR-T model.
#'
#' @param y Array with the time series observations
#' @param sT Array with the irregular observational times
#' @param nu degrees of freedom
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR-T model.}
#' \item{sigma}{ MLE of the sigma parameter of the IAR-T model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi and sigma.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARt.sample}}, \code{\link{IAR.phi.t}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-IARt.sample(n,0.9,st,sigma2=1,nu=3)
#' model<-IAR.t(y$y, sT=st)
#' phi=model$phi
#' sigmaest=model$sigma
IAR.t<-function (y, sT,nu=3) #Find minimum of IAR.phi.gamma
{
  aux<-1e10
  value<-1e10
  br<-0
  for(i in 1:20)
  {
    phi=runif(1)
    sigma=var(y)*runif(1)
    optim<-nlminb(start=c(phi,sigma),objective=IAR.phi.t,y=y,sT=sT,nu=nu,lower=c(0,0.0001),upper=c(0.9999,2*var(y)))
    value<-optim$objective
    #print(c(optim$objective,optim$par,aux>value))
    if(aux>value)
    {
      par<-optim$par
      aux<-value
      br<-br+1
    }
    if(aux<=value & br>10 & i>15)
      break;
  }
  if(aux==1e10)
    par<-c(0,0)
  return(list(phi=par[1],sigma=par[2],ll=aux))
}
