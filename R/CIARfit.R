#' Fitted Values of CIAR model
#'
#' Fit a CIAR model to an irregularly observed time series.
#'
#' @param x An array with the parameters of the CIAR model. The elements of the array are, in order, the real and the imaginary part of the phi parameter of the CIAR model.
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{yhat}{ Fitted values of the observable part of CIAR model.}
#' \item{xhat}{ Fitted values of both observable part and imaginary part of CIAR model.}
#' \item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
#' \item{Theta}{ Theta array estimated by the CIAR model at the last time point.}
#' \item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
#' \item{Qt}{ Covariance matrix of the state equation estimated by the CIAR model at the last time point.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIAR.sample}}, \code{\link{CIAR.phi.kalman}},\code{\link{CIAR.kalman}}
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
#' yhat=CIAR.fit(x=c(ciar$phiR,ciar$phiI),y=y1,t=st)
CIAR.fit<-function(x,y,t,standarized='TRUE',c=1)
{
  sigmay<-1
  if(standarized=='FALSE')
    sigmay<-var(y)
  n=length(y)
  Sighat=sigmay*matrix(c(1,0,0,c),2,2)
  xhat=matrix(0,nrow=2,ncol=n)
  delta<-diff(t)
  Q=Sighat
  phi.R=x[1]
  phi.I=x[2]
  #DEFINITION OF F
  F=matrix(0,nrow=2,ncol=2)
  G=matrix(0,nrow=1,ncol=2)
  G[1,1]=1
  #MOD PHI MUST BE LESS THAN ONE
  phi=complex(1,real=phi.R,imaginary=phi.I)
  Phi=Mod(phi)
  phi=ifelse(is.na(phi)==TRUE,1.1,phi)
  if(Mod(phi)>=1)
    stop("Mod of Phi must be less than one")
  psi<-acos(phi.R/Phi)
  if(Mod(phi)<1){
    for(i in 1:(n-1))
    {
      Lambda=G%*%Sighat%*%t(G) #R_t es 0
      if(Lambda<=0 | is.na(Lambda)==TRUE)
        break;
      phi2.R<-(Phi**delta[i])*cos(delta[i]*psi)
      phi2.I<-(Phi**delta[i])*sin(delta[i]*psi)
      F[1,1]=phi2.R
      F[1,2]=-phi2.I
      F[2,1]=phi2.I
      F[2,2]=phi2.R
      phi2<-1-Mod(phi**delta[i])**2
      Qt<-phi2*Q
      Theta=F%*%Sighat%*%t(G)
      xhat[,i+1]=F%*%xhat[,i]+Theta%*%solve(Lambda)%*%(y[i]-G%*%xhat[,i])
      Sighat=F%*%Sighat%*%t(F)+ Qt - Theta%*%solve(Lambda)%*%t(Theta)
    }
    yhat=G%*%xhat
  }
  return(list(yhat=yhat,xhat=xhat,Sighat=Sighat,Theta=Theta,Lambda=Lambda,Qt=Qt))
}
