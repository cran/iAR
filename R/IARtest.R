#' Test for the significance of the autocorrelation estimated by the IAR model in periodic irregularly observed time series
#'
#' This function perform a test for the significance of the autocorrelation estimated by the IAR model. This test is based on the residuals of the periodical time series fitted with an harmonic model using an incorrect period.
#'
#' @param y Array with the time series observations
#' @param st Array with the irregular observational times
#' @param merr Array with the variance of the measurement errors.
#' @param f Frequency (1/Period) of the raw time series
#' @param phi autocorrelation estimated by \code{\link{IARloglik}}
#' @param plot logical; if true, the function return a density plot of the distribution of the bad fitted examples; if false, this function does not return a plot
#' @param xlim The x-axis limits (x1, x2) of the plot. Only works if plot='TRUE'. See \code{\link{plot.default}} for more details
#'
#' @details The null hypothesis of the test is: The autocorrelation estimated in the time series belongs to the distribution of the coefficients estimated for the residuals of the data fitted using wrong periods. Therefore, if the hypothesis is rejected, it can be concluded that the residuals of the harmonic model do not remain a time dependency structure.The statistic of the test is log(phi) which was contrasted with a normal distribution with parameters corresponding to the log of the mean and the variance of the phi computed for the residuals of the bad fitted light curves.
#' @return A list with the following components:
#' \itemize{
#' \item{phi}{ MLE of the phi parameter of the IAR model.}
#' \item{bad}{ A matrix with two columns. The first column contains the incorrect frequencies used to fit each harmonic model. The second column has the MLEs of the phi parameters of the IAR model that has been fitted to the residuals of the harmonic model fitted using the frequencies of the first column.}
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
#' \code{\link{clcep}}, \code{\link{harmonicfit}}, \code{\link{IARloglik}}, \code{\link{IARTest2}}
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' results=harmonicfit(file=clcep,f1=f1)
#' y=results$res/sqrt(var(results$res))
#' st=results$t
#' res3=IARloglik(y,st,standarized='TRUE')[1]
#' res3$phi
#' require(ggplot2)
#' test<-IARTest(y=clcep[,2],st=clcep[,1],f=f1,phi=res3$phi,plot='TRUE',xlim=c(-10,0.5))
#' test
#' outbad <- ggplot(test$bad, aes(x = f0, y = bad)) + geom_line()+geom_point() +
#' geom_point(aes(x=f1,y=as.numeric(test$phi)),color="red") +
#' xlab("Frequencies") + ylab("iAR Coefficient") +
#' ylim(c(0,max(c(as.numeric(test$phi),test$bad$bad))))+theme_bw() +
#' ggtitle("")+
#' theme(plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
#' panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#' panel.background = element_blank())
#' outbad
IARTest<-function (y, st, merr=0,f, phi, plot = T, xlim = c(-1, 0))
{
  aux = seq(2.5, 47.5, by = 2.5)
  aux = c(-aux, aux)
  aux = sort(aux)
  f0 = f * (1 + aux/100)
  f0 <- sort(f0)
  l1 <- length(f0)
  bad <- rep(0, l1)
  data <- cbind(st, y)
  if(sum(merr)==0)
    merr=rep(0,length(y))
  for (j in 1:l1) {
    results = harmonicfit(file = data, f1 = f0[j])
    y = results$res/sqrt(var(results$res))
    st = results$t
    res3 = IARloglik(y, st,merr)[1]
    bad[j] = res3$phi
  }
  mubf <- mean(log(bad))
  sdbf <- sd(log(bad))
  z0 <- log(phi)
  pvalue <- pnorm(z0, mubf, sdbf)
  out <- NULL
  BAD <-data.frame(f0,bad)
  if (plot == T) {
    phi2 = bad
    phi2 <- as.data.frame(phi2)
    phi <- as.data.frame(phi)
    out <- ggplot(phi2, aes(log(phi2))) + geom_density(adjust = 2) +
      xlab("") + ylab("") + theme_bw() + ggtitle("") +
      geom_point(data = phi, aes(log(phi)), y = 0, size = 4,
                 col = "red", shape = 17) + xlim(xlim[1], xlim[2]) +
      theme(plot.title = element_text(face = "bold", size = 20),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    out
  }
  return(list(phi = phi, bad=BAD, norm = c(mubf, sdbf), z0 = z0, pvalue = pvalue, out = out))
}
