#' Maximum Likelihood Estimation of the IAR Model
#'
#' Maximum Likelihood Estimation of the IAR Model.
#'
#' @param y Array with the time series observations.
#' @param st Array with the irregular observational times.
#' @param delta Array with the measurements error standard deviations.
#' @param zero.mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi.}
#' }
#' @export
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{IARsample}}, \code{\link{arima}}, \code{\link{IARphiloglik}}
#'
#' @examples
#' #Generating IAR sample
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-IARsample(phi=0.99,st=st,n=100)
#' y<-y$series
#' #Compute Phi
#' phi=IARloglik(y=y,st=st)$phi
#' print(phi)
#' #Compute the standard deviation of innovations
#' n=length(y)
#' d=c(0,diff(st))
#' phi1=phi**d
#' yhat=phi1*as.vector(c(0,y[1:(n-1)]))
#' plot(st,y,type='l')
#' lines(st,yhat,col='red')
#' sigma=var(y)
#' nu=c(sigma,sigma*(1-phi1**(2))[-1])
#' tau<-nu/sigma
#' sigmahat<-mean(c((y-yhat)**2/tau))
#' nuhat<-sigmahat*(1-phi1**(2))
#' nuhat2<-sqrt(nuhat)
#' #Equally spaced models
#' require(arfima)
#' fit2<-arfima(y,order=c(1,0,0))
#' fit<-arima(y,order=c(1,0,0),include.mean=FALSE)
#' syarf<-tacvfARFIMA(phi=fit2$modes[[1]]$phi,dfrac=fit2$modes[[1]]$dfrac,
#' sigma2=fit2$modes[[1]]$sigma,maxlag=20)[1]
#' syar<-fit$sigma/(1-fit$coef[1]**2)
#' print(sigmahat)
#' print(syar)
#' print(syarf)
#' carf<-fit2$modes[[1]]$sigma/syarf
#' car<-(1-fit$coef[1]**2)
#' ciar<-(1-phi1**(2))
#' #Compute the standard deviation of innovations (regular case)
#' sigma=var(y)
#' nuhat3=sqrt(sigma*ciar)
#' searf<-sqrt(sigma*carf)
#' sear<-sqrt(sigma*car)
#' #Plot the standard deviation of innovations
#' plot(st[-1], nuhat3[-1], t="n", axes=FALSE,xlab='Time',ylab='Standard Deviation of Innovations')
#' axis(1)
#' axis(2)
#' segments(x0=st[-1], y0=nuhat3[-1], y1=0, col=8)
#' points(st, nuhat3, pch=20, col=1, bg=1)
#' abline(h=sd(y),col='red',lwd=2)
#' abline(h=sear,col='blue',lwd=2)
#' abline(h=searf,col='green',lwd=2)
#' abline(h=mean(nuhat3[-1]),col='black',lwd=2)
IARloglik=function(y,st,delta=0,zero.mean=TRUE,standardized=TRUE){
  if(sum(delta)==0){
    delta=rep(0,length(y))}
  out=optimize(IARphiloglik, interval=c(0,1), y=y, st=st, delta_input=delta, zeroMean = zero.mean, standardized = standardized)
  phi=out$minimum
  ll=out$objective
  return(list(phi=phi,loglik=ll))
}
