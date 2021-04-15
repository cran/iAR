#' Minus Log Likelihood of the CIAR Model
#'
#' This function return the negative log likelihood of the CIAR process given specific values of phi.R and phi.I
#'
#' @param x An array with the parameters of the CIAR model. The elements of the array are, in order, the real (phi.R) and the imaginary (phi.I) part of the coefficient of CIAR model.
#' @param y Array with the time series observations.
#' @param t Array with the irregular observational times.
#' @param yerr Array with the measurements error standard deviations.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @return Value of the negative log likelihood evaluated in phiR and phiI.
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CIAR.sample}}
#'
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIAR.sample(n=n,phi.R=0.9,phi.I=0,sT=st,c=1)
#' y=x$y
#' yerr=rep(0,n)
#' CIAR.phi.kalman(x=c(0.8,0),y=y,t=st,yerr=yerr)
CIAR.phi.kalman<-function(x,y,t,yerr,zero.mean='TRUE',standarized='TRUE',c=1)
{
  sigmay<-1
  if(standarized=='FALSE')
    sigmay<-var(y)
  if(zero.mean=='FALSE')
    y=y-mean(y)
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
  psi<-acos(phi.R/Phi)
  sum.Lambda=0
  sum.error=0
  phi=ifelse(is.na(phi)==TRUE,1.1,phi)
  if(Mod(phi)<1){
    for(i in 1:(n-1))
    {
      Lambda=G%*%Sighat%*%t(G) + yerr[i+1]**2
      if(Lambda<=0 | is.na(Lambda)==TRUE)
      {
        sum.Lambda<-n*1e10
        break;
      }
      phi2.R<-(Phi**delta[i])*cos(delta[i]*psi)
      phi2.I<-(Phi**delta[i])*sin(delta[i]*psi)
      F[1,1]=phi2.R
      F[1,2]=-phi2.I
      F[2,1]=phi2.I
      F[2,2]=phi2.R
      phi2<-1-Mod(phi**delta[i])**2
      Qt<-phi2*Q
      sum.Lambda=sum.Lambda+log(Lambda)
      Theta=F%*%Sighat%*%t(G)
      sum.error= sum.error+ ( y[i]-G%*%xhat[,i] )**2/Lambda
      xhat[,i+1]=F%*%xhat[,i]+Theta%*%solve(Lambda)%*%(y[i]-G%*%xhat[,i])
      Sighat=F%*%Sighat%*%t(F)+ Qt - Theta%*%solve(Lambda)%*%t(Theta)
    }
    yhat=G%*%xhat
    out<-ifelse(is.na(sum.Lambda)==TRUE,1e10,(sum.Lambda + sum.error)/n)
  }
  else out=1e10
  return(out)
}
