#' Minus Log Likelihood of the IAR Model estimated via Kalman Recursions
#'
#' This function return the negative log likelihood of the IAR process given a specific value of phi.
#'
#' @param x A given phi coefficient of the IAR model.
#' @param y Array with the time series observations.
#' @param yerr Array with the measurements error standard deviations.
#' @param t Array with the irregular observational times.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
#'
#' @return Value of the negative log likelihood evaluated in phi.
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IAR.sample}}
#'
#' @examples
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IAR.sample(phi=0.99,n=100,st)
#' y<-y$series
#' IAR.phi.loglik(x=0.8,y=y,sT=st)
IAR.phi.kalman<-function(x,y,yerr,t,zero.mean='TRUE',standarized='TRUE')
{
  sigmay<-1
  if(standarized=='FALSE')
    sigmay<-var(y)
  if(zero.mean=='FALSE')
    y=y-mean(y)
  n=length(y)
  Sighat=sigmay*diag(1)
  xhat=matrix(0,nrow=1,ncol=n)
  delta<-diff(t)
  Q=Sighat
  phi=x
  #DEFINITION OF F
  F=matrix(0,nrow=1,ncol=1)
  G=matrix(0,nrow=1,ncol=1)
  G[1,1]=1
  #MOD PHI MUST BE LESS THAN ONE
  sum.Lambda=0
  sum.error=0
  phi=ifelse(is.na(phi)==TRUE,1.1,phi)
  if(abs(phi)<1){
    for(i in 1:(n-1))
    {
      Lambda=G%*%Sighat%*%t(G) + yerr[i+1]**2
      if(Lambda<=0 | is.na(Lambda)==TRUE)
      {
        sum.Lambda<-n*1e10
        break;
      }
      phi2<-phi**delta[i]
      F[1,1]=phi2
      phi2<-1-phi**(2*delta[i])
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
