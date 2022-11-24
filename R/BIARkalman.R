#' Maximum Likelihood Estimation of the BIAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the BIAR model parameters phiR and phiI. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#'
#' @param y1 Array with the observations of the first time series of the BIAR process.
#' @param y2 Array with the observations of the second time series of the BIAR process.
#' @param t Array with the irregular observational times.
#' @param delta1 Array with the measurements error standard deviations of the first time series of the BIAR process.
#' @param delta2 Array with the measurements error standard deviations of the second time series of the BIAR process.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param niter Number of iterations in which the function nlminb will be repeated.
#' @param seed a single value, interpreted as the seed of the random process.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phiR}{ MLE of the autocorrelation coefficient of BIAR model (phiR).}
#' \item{phiI}{ MLE of the cross-correlation coefficient of the BIAR model (phiI).}
#' \item{ll}{ Value of the negative log likelihood evaluated in phiR and phiI.}
#' }
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{BIARsample}}, \code{\link{BIARphikalman}}
#'
#'
#' @examples
#' \donttest{
#' n=80
#' set.seed(6714)
#' st<-gentime(n)
#' x=BIARsample(n=n,phiR=0.9,phiI=0,st=st,rho=0)
#' y=x$y
#' y1=y/apply(y,1,sd)
#' biar=BIARkalman(y1=y1[1,],y2=y1[2,],t=st,delta1 = rep(0,length(y[1,])),
#' delta2=rep(0,length(y[1,])))
#' biar
#' }
BIARkalman<-function (y1, y2, t, delta1 = 0, delta2 = 0, zero.mean = "TRUE", niter = 10, seed = 1234)
{
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
  for (i in 1:niter) {
    phiR = 2 * runif(1) - 1
    phiI = 2 * runif(1) - 1
    if (Mod(complex(1, real = phiR, imaginary = phiI)) < 1) {
      optim <- nlminb(start = c(phiR, phiI), objective = BIARphikalman,
                      y1 = y1,y2 = y2, t = t, yerr1 = delta1, yerr2=delta2,yest=c(0,0),lower	= c(-1, -1), upper = c(1, 1))
      value <- optim$objective
    }
    if (aux > value) {
      par <- optim$par
      aux <- value
      br <- br + 1
    }
    if (aux <= value & br > 1 & i > trunc(niter/2))
      break
  }
  if (aux == 1e+10)
    par <- c(0, 0)
  return(list(phiR = par[1], phiI = par[2], ll = aux))
}
