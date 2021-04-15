#' Test for the significance of the autocorrelation estimated by the IAR model
#'
#' This function perform a test for the significance of the autocorrelation estimated by the IAR model. This test is based in to take N disordered samples of the original data (Useful for non-periodic time series or when the period is unknown).
#'
#' @param y Array with the time series observations
#' @param sT Array with the irregular observational times
#' @param iter Number of disordered samples of the original data (N)
#' @param phi autocorrelation estimated by \code{\link{IAR.loglik}}
#' @param plot logical; if true, the function return a density plot of the distribution of the bad fitted examples; if false, this function does not return a plot
#' @param xlim The x-axis limits (x1, x2) of the plot. Only works if plot='TRUE'. See \code{\link{plot.default}} for more details
#'
#' @details The null hypothesis of the test is: The coefficient phi estimated for the time series belongs to the distribution of the coefficients estimated on the disordered data, which are assumed to be uncorrelated. Therefore, if the hypothesis is accepted, it can be concluded that the observations of the time series are uncorrelated.The statistic of the test is log(phi) which was contrasted with a normal distribution with parameters corresponding to the log of the mean and the variance of the phi computed for the N samples of the disordered data. This test differs for \code{\link{IAR.Test}} in that to perform this test it is not necessary to know the period of the time series.
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR model.}
#' \item{norm}{ Mean and variance of the normal distribution of the disordered data.}
#' \item{z0}{ Statistic of the test (log(phi)).}
#' \item{pvalue}{ P-value computed for the test.}
#' }
#' @export
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{Planets}},\code{\link{IAR.loglik}}, \code{\link{IAR.Test}}
#'
#' @examples
#' data(Planets)
#' t<-Planets[,1]
#' res<-Planets[,2]
#' y=res/sqrt(var(res))
#' res3=IAR.loglik(y,t,standarized='TRUE')[1]
#' res3$phi
#' set.seed(6713)
#' require(ggplot2)
#' test<-IAR.Test2(y=y,sT=t,phi=res3$phi,plot='TRUE',xlim=c(-9.6,-9.45))
IAR.Test2=function(y,sT,iter=100,phi,plot='TRUE',xlim=c(-1,0))
{
  phi2=rep(0,iter)
  for(i in 1:iter)
  {
    ord<-sample(1:length(y))
    y1<-y[ord]
    phi2[i]=IAR.loglik(y=y1,sT=sT)$phi
  }
  mubf<-mean(log(phi2))
  sdbf<-sd(log(phi2))
  z0<-log(phi)
  pvalue<-pnorm(z0,mubf,sdbf)
  out<-NULL
  if(plot=='TRUE')
  {
    phi2<-as.data.frame(phi2)
    phi<-as.data.frame(phi)
    out<-ggplot(phi2,aes(log(phi2)))+geom_density(adjust=2)+xlab("")		+ylab("")+theme_bw()+ggtitle("")+geom_point(data = phi,aes(log(phi)), y = 0, size = 4,col='red',shape=17) + xlim(xlim[1],xlim[2])+
      theme(plot.title = element_text(face="bold", size=20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    out
  }
  return(list(phi=phi,norm=c(mubf,sdbf),z0=z0,pvalue=pvalue,out=out))
}
