#' Smoothing and Forecasting from BIAR model
#'
#' Estimation of missing values from models fitted by \code{\link{BIARkalman}}
#'
#' @param x An array with the parameters of the BIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of BIAR model.
#' @param y1 Array with the observations of the first time series of the BIAR process.
#' @param y2 Array with the observations of the second time series of the BIAR process.
#' @param t Array with the irregular observational times.
#' @param delta1 Array with the measurements error standard deviations of the first time series of the BIAR process.
#' @param delta2 Array with the measurements error standard deviations of the second time series of the BIAR process.
#' @param yini1 a single value, initial value of the estimation of the missing value of the first time series of the BIAR process.
#' @param yini2 a single value, initial value of the estimation of the missing value of the second time series of the BIAR process.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param niter Number of iterations in which the function nlminb will be repeated.
#' @param seed a single value, interpreted as the seed of the random process.
#' @param nsmooth a single value; If 1, only one time series of the BIAR process has a missing value. If 2, both time series of the BIAR process have a missing value.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{Estimation of the missing values of the BIAR process.}
#' \item{ll}{Value of the negative log likelihood evaluated in the fitted missing values.}
#' }

#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{BIARsample}}, \code{\link{BIARLL}}
#'
#' @examples
#' \donttest{
#' set.seed(6713)
#' n=100
#' st<-gentime(n)
#' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st,rho=0.9)
#' y=x$y
#' y1=y/apply(y,1,sd)
#' yerr1=rep(0,n)
#' yerr2=rep(0,n)
#' biar=BIARkalman(y1=y1[1,],y2=y1[2,],t=st,delta1 = yerr1,delta2=yerr2)
#' biar
#' napos=10
#' y0=y1
#' y1[1,napos]=NA
#' xest=c(biar$phiR,biar$phiI)
#' yest=BIARsmoothing(xest,y1=y1[1,],y2=y1[2,],t=st,delta1=yerr1,
#' delta2=yerr2,nsmooth=1)
#' yest$fitted
#' mse=(y0[1,napos]-yest$fitted)^2
#' print(mse)
#' par(mfrow=c(2,1))
#' plot(st,x$y[1,],type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,x$y[1,],pch=20)
#' points(st[napos],yest$fitted*apply(y,1,sd)[1],col="red",pch=20)
#' plot(st,x$y[2,],type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,x$y[2,],pch=20)
#' }
BIARsmoothing <- function (x,y1, y2, t, delta1 = 0, delta2 = 0, yini1=0,yini2=0,zero.mean = "TRUE", niter = 10, seed = 1234,nsmooth=1) {
  set.seed(seed)
  aux <- 1e+10
  value <- 1e+10
  br <- 0

  if (sum(delta1) == 0) {
    delta1 = rep(0, length(y1))
  }

  if (sum(delta2) == 0) {
    delta2 = rep(0, length(y2))
  }

  if (yini1==0)
    yini1 = rnorm(1)

  if (yini2==0)
    yini2 = rnorm(1)

  if (nsmooth==2){
    for (i in 1:niter) {
      optim <- nlminb(start = c(yini1, yini2), objective = BIARLL,
                      phiValues=x,y1 = y1,y2 = y2, t = t, yerr1 = delta1,
                      yerr2=delta2, zeroMean = zero.mean,
                      lower = c(-Inf, -Inf), upper = c(Inf, Inf))

      value <- optim$objective

      if (aux > value) {
        par <- optim$par
        aux <- value
        br <- br + 1
      }

      if (aux <= value & br > 1 & i > trunc(niter/2))
        break
    }
  }

  if (aux == 1e+10)
    par <- c(0, 0)

  if (nsmooth==1){
    out = optim(c(yini1), BIARLL, lower = -Inf, upper = Inf,
                phiValues=x, y1=y1, y2=y2, t=t, yerr1 = delta1,
                yerr2=delta2, zeroMean = zero.mean, method="BFGS")
    par = out$par
    aux = out$value
  }

  return(list(fitted = par, ll = aux))
}
