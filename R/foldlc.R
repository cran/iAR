#' Plotting folded light curves
#'
#' This function plots a time series folded on its period.
#' @param file Matrix with the light curve observations. The first column must have the irregular times, the second column must have the brightness magnitudes and the third column must have the measurement errors.
#' @param f1 Frequency (1/Period) of the light curve.
#'
#' @return A plot of the folded (phased) time series.
#' @export
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' foldlc(clcep,f1)
foldlc<-function(file,f1)
{
  mycurve=file
  t=mycurve[,1]
  m=mycurve[,2]
  merr=mycurve[,3]
  P<-1/f1
  fold<-(t-t[1])%%(P)/P
  fold<-c(fold,fold+1)
  m<-rep(m,2)
  merr<-rep(merr,2)
  dat1<-cbind(fold,m,merr)
  dat1<-as.data.frame(dat1)
  out<-ggplot(dat1, aes(x=fold, y=m)) +
    geom_errorbar(aes(ymin=m-merr, ymax=m+merr), width=.01,col='red') +
    geom_point()+scale_y_reverse()+xlab("")+ylab("")+theme_bw()+
    theme(plot.title = element_text(face="bold", size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  out
}
