#' Test for the significance of the autocorrelation estimated by the IAR model in periodic irregularly observed time series
#'
#' This function perform a test for the significance of the autocorrelation estimated by the IAR model. This test is based on the residuals of the periodical time series fitted with an harmonic model using an incorrect period.
#'
#' @param y Array with the time series observations
#' @param sT Array with the irregular observational times
#' @param f Frequency (1/Period) of the raw time series
#' @param phi autocorrelation estimated by \code{\link{IAR.loglik}}
#' @param plot logical; if true, the function return a density plot of the distribution of the bad fitted examples; if false, this function does not return a plot
#' @param xlim The x-axis limits (x1, x2) of the plot. Only works if plot='TRUE'. See \code{\link{plot.default}} for more details
#'
#' @details The null hypothesis of the test is: The autocorrelation estimated in the time series belongs to the distribution of the coefficients estimated for the residuals of the data fitted using wrong periods. Therefore, if the hypothesis is rejected, it can be concluded that the residuals of the harmonic model do not remain a time dependency structure.The statistic of the test is log(phi) which was contrasted with a normal distribution with parameters corresponding to the log of the mean and the variance of the phi computed for the residuals of the bad fitted light curves.
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR model.}
#' \item{norm}{ Mean and variance of the normal distribution of the bad fitted examples.}
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
#' \code{\link{clcep}}, \code{\link{harmonicfit}}, \code{\link{IAR.loglik}}, \code{\link{IAR.Test2}}
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' results=harmonicfit(file=clcep,f1=f1)
#' y=results$res/sqrt(var(results$res))
#' sT=results$t
#' res3=IAR.loglik(y,sT,standarized='TRUE')[1]
#' res3$phi
#' require(ggplot2)
#' test<-IAR.Test(y=clcep[,2],sT=clcep[,1],f1,res3$phi,plot='TRUE',xlim=c(-10,0.5))
#' test
IAR.Test=function(y,sT,f,phi,plot='TRUE',xlim=c(-1,0))
{
  aux=seq(2.5,47.5,by=2.5)
  aux=c(-aux,aux)
  aux=sort(aux)
  f0=f*(1+aux/100)
  f0<-sort(f0)
  l1<-length(f0)
  bad<-rep(0,l1)
  data<-cbind(sT,y)
  for(j in 1:l1)
  {
    results=harmonicfit(file=data,f1=f0[j])
    y=results$res/sqrt(var(results$res))
    sT=results$t
    res3=IAR.loglik(y,sT)[1]
    bad[j]=res3$phi
  }
  mubf<-mean(log(bad))
  sdbf<-sd(log(bad))
  z0<-log(phi)
  pvalue<-pnorm(z0,mubf,sdbf)
  out<-NULL
  if(plot=='TRUE')
  {
    phi2=bad
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
