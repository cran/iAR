#' Simulate from a CIAR Model
#'
#' Simulates a CIAR Time Series Model
#'
#' @param n Length of the output time series. A strictly positive integer.
#' @param sT Array with observational times.
#' @param phi.R Real part of the coefficient of CIAR model. A value between -1 and 1.
#' @param phi.I Imaginary part of the coefficient of CIAR model. A value between -1 and 1.
#' @param rho Correlation between the real and the imaginary part of the process. A value between -1 and 1.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @details The chosen phi.R and phi.I values must satisfy the condition $|phi.R + i phi.I| < 1$.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{y}{Array with the simulated real part of the CIAR process.}
#' \item{t}{ Array with observation times.}
#' \item{Sigma}{ Covariance matrix of the process.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIAR.sample(n=n,phi.R=0.9,phi.I=0,sT=st,c=1)
#' plot(st,x$y,type='l')
#' x=CIAR.sample(n=n,phi.R=-0.9,phi.I=0,sT=st,c=1)
#' plot(st,x$y,type='l')
CIAR.sample<-function(n,sT,phi.R,phi.I,rho=0,c=1)
{
  delta<-diff(sT)
  x=matrix(0,nrow=2,ncol=n)
  F=matrix(0,nrow=2,ncol=2)
  phi=complex(1,real=phi.R,imaginary=phi.I)
  if(Mod(phi)>=1)
    stop("Mod of Phi must be less than one")
  Phi=Mod(phi)
  psi<-acos(phi.R/Phi)
  e.R=rnorm(n)
  e.I=rnorm(n)
  state.error=rbind(e.R,e.I)
  Sigma=matrix(1,nrow=2,ncol=2)
  Sigma[1,1]=1
  Sigma[2,2]=c
  Sigma[1,2]=rho*sqrt(Sigma[1,1])*sqrt(Sigma[2,2])
  Sigma[2,1]=Sigma[1,2]
  B=svd(Sigma)
  A=matrix(0,nrow=2,ncol=2)
  diag(A)=sqrt(B$d)
  Sigma.root=(B$u)%*%A%*%B$u
  state.error=Sigma.root%*%state.error
  G=matrix(0,nrow=1,ncol=2)
  G[1,1]=1
  y=numeric()
  x[,1]=state.error[,1]
  for(i in 1:(n-1))
  {
    phi2.R<-(Phi**delta[i])*cos(delta[i]*psi)
    phi2.I<-(Phi**delta[i])*sin(delta[i]*psi)
    phi2<-1-Mod(phi**delta[i])**2
    F[1,1]=phi2.R
    F[1,2]=-phi2.I
    F[2,1]=phi2.I
    F[2,2]=phi2.R
    x[,i+1]=F%*%x[,i]+sqrt(phi2)*state.error[,i]
    y[i]=G%*%x[,i]
  }
  y[n]=G%*%x[,n]
  return(list(t=sT,y=y,Sigma=Sigma))
}
