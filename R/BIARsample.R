#' Simulate from a BIAR Model
#'
#' Simulates a BIAR Time Series Model
#'
#' @param n Length of the output bivariate time series. A strictly positive integer.
#' @param st Array with observational times.
#' @param phiR Autocorrelation coefficient of BIAR model. A value between -1 and 1.
#' @param phiI Crosscorrelation coefficient of BIAR model. A value between -1 and 1.
#' @param delta1 Array with the measurements error standard deviations of the first time series of the bivariate process.
#' @param delta2 Array with the measurements error standard deviations of the second time series of the bivariate process.
#' @param rho Contemporary correlation coefficient of BIAR model. A value between -1 and 1.
#'
#' @details The chosen phiR and phiI values must satisfy the condition $|phiR + i phiI| < 1$.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{y}{ Matrix with the simulated BIAR process.}
#' \item{t}{ Array with observation times.}
#' \item{Sigma}{ Covariance matrix of the process.}
#' }
#'
#' @export
#'
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#' @seealso
#' \code{\link{gentime}}
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
#' plot(st,x$y[1,],type='l')
#' plot(st,x$y[2,],type='l')
#' x=BIARsample(n=n,phiR=-0.9,phiI=-0.3,st=st)
#' plot(st,x$y[1,],type='l')
#' plot(st,x$y[2,],type='l')
BIARsample<-function (n, st, phiR, phiI, delta1=0,delta2=0,rho = 0)
{
  delta <- diff(st)
  x = matrix(0, nrow = 2, ncol = n)
  F = matrix(0, nrow = 2, ncol = 2)
  phi = complex(1, real = phiR, imaginary = phiI)
  if (Mod(phi) >= 1)
    stop("Mod of Phi must be less than one")
  Phi = Mod(phi)
  psi <- acos(phiR/Phi)
  if (phiI < 0)
    psi = -acos(phiR/Phi)
  #psi2 <- asin(phiI/Phi)
  e.R = rnorm(n)
  e.I = rnorm(n)
  state.error = rbind(e.R, e.I)
  Sigma = matrix(1, nrow = 2, ncol = 2)
  Sigma[1, 1] = 1
  Sigma[2, 2] = 1
  Sigma[1, 2] = rho * sqrt(Sigma[1, 1]) * sqrt(Sigma[2, 2])
  Sigma[2, 1] = Sigma[1, 2]
  B = svd(Sigma)
  A = matrix(0, nrow = 2, ncol = 2)
  diag(A) = sqrt(B$d)
  Sigma.root = (B$u) %*% A %*% B$u
  state.error = Sigma.root %*% state.error
  #measurement error
  w.R = rnorm(n)
  w.I = rnorm(n)
  observation.error = rbind(w.R, w.I)
  SigmaO = matrix(0, nrow = 2, ncol = 2)
  SigmaO[1, 1] = delta1**2
  SigmaO[2, 2] = delta2**2
  BO = svd(SigmaO)
  AO = matrix(0, nrow = 2, ncol = 2)
  diag(AO) = sqrt(BO$d)
  SigmaO.root = (BO$u) %*% AO %*% BO$u
  observation.error = SigmaO.root %*% observation.error
  G = diag(2)
  y = observation.error
  x[, 1] = state.error[, 1]
  for (i in 1:(n - 1)) {
    phi2.R <- (Phi^delta[i]) * cos(delta[i] * psi)
    phi2.I <- (Phi^delta[i]) * sin(delta[i] * psi)
    phi2 <- 1 - Mod(phi^delta[i])^2
    F[1, 1] = phi2.R
    F[1, 2] = -phi2.I
    F[2, 1] = phi2.I
    F[2, 2] = phi2.R
    x[, i + 1] = F %*% x[, i] + sqrt(phi2) * state.error[,i]
    y[, i ] = G %*% x[, i] + observation.error[,i]
  }
  y[, n ] = G %*% x[, n] + observation.error[,n]
  return(list(t = st, y = y, Sigma = Sigma))
}
