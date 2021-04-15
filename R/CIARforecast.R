#' Forecast from CIAR model
#'
#' Forecast from models fitted by \code{\link{CIAR.kalman}}
#'
#' @param phi.R Real part of the phi coefficient of CIAR model.
#' @param phi.I Imaginary part of the phi coefficient of CIAR model.
#' @param y1 Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param n.ahead The number of steps ahead for forecast is required.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Fitted values by the CIAR model.}
#' \item{forecast}{ Point Forecasts in the n.ahead times.}
#' \item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
#' \item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{CIAR.sample}}, \code{\link{CIAR.kalman}}, \code{\link{CIAR.fit}}
#'
#'
#' @examples
#' #Simulated Data
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=CIAR.sample(n=n,phi.R=0.9,phi.I=0,sT=st,c=1)
#' y=x$y
#' y1=y/sd(y)
#' n=length(y1)
#' p=trunc(n*0.99)
#' ytr=y1[1:p]
#' yte=y1[(p+1):n]
#' str=st[1:p]
#' ste=st[(p+1):n]
#' n.ahead=ste-str[p]
#'
#' final<-matrix(0,length(n.ahead),4)
#' ciar=CIAR.kalman(y=ytr,t=str)
#' forCIAR<-CIAR.forecast(ciar$phiR,ciar$phiI,ytr,str,n.ahead=n.ahead)
#'
CIAR.forecast<-function(phi.R,phi.I,y1,st,n.ahead=1)
{
  yhat=CIAR.fit(x=c(phi.R,phi.I),y=y1,t=st)
  n=length(yhat$yhat)
  xhat=yhat$xhat
  Theta=yhat$Theta
  Lambda=yhat$Lambda
  Sighat=yhat$Sighat
  Qt=yhat$Qt
  F=matrix(0,nrow=2,ncol=2)
  G=matrix(0,nrow=1,ncol=2)
  G[1,1]=1
  phi=complex(1,real=phi.R,imaginary=phi.I)
  Phi=Mod(phi)
  phi=ifelse(is.na(phi)==TRUE,1.1,phi)
  if(Mod(phi)>=1)
    stop("Mod of Phi must be less than one")
  psi<-acos(phi.R/Phi)
  delta=n.ahead
  phi2.R<-(Phi**delta)*cos(delta*psi)
  phi2.I<-(Phi**delta)*sin(delta*psi)
  yhat1=rep(0,length(n.ahead))
  xhat1=matrix(0,nrow=2,ncol=length(n.ahead))
  Lambda2=rep(0,length(n.ahead))
  for(i in 1:length(n.ahead))
  {
    F[1,1]=phi2.R[i]
    F[1,2]=-phi2.I[i]
    F[2,1]=phi2.I[i]
    F[2,2]=phi2.R[i]
    xhat1[,i]=F%*%xhat[,n] +Theta%*%solve(Lambda)%*%(y1[n]-G%*%xhat[,n])
    yhat1[i]=G%*%xhat1[,i]
    Sighat2=F%*%Sighat%*%t(F)+ Qt - Theta%*%solve(Lambda)%*%t(Theta) #R_t es 0
    Lambda2[i]=G%*%Sighat2%*%t(G)
  }
  return(list(fitted=yhat$yhat,forecast=yhat1,Lambda=Lambda2,Sighat=Sighat2))
}
