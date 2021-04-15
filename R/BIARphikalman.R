#' Minus Log Likelihood of the BIAR Model
#'
#' This function return the negative log likelihood of the BIAR process given specific values of phi.R and phi.I
#'
#'
#' @param x An array with the parameters of the BIAR model. The elements of the array are, in order, the real (phi.R) and the imaginary (phi.I) part of the coefficient of BIAR model.
#' @param y1 Array with the observations of the first time series of the BIAR process.
#' @param y2 Array with the observations of the second time series of the BIAR process.
#' @param t Array with the irregular observational times.
#' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
#' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
#' @param zero.mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#'
#' @return Value of the negative log likelihood evaluated in phiR and phiI.
#' @export
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{BIAR.sample}}
#'
#'
#' @examples
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' x=BIAR.sample(n=n,phi.R=0.9,phi.I=0.3,sT=st)
#' y=x$y
#' y1=y[1,]
#' y2=y[2,]
#' yerr1=rep(0,n)
#' yerr2=rep(0,n)
#' BIAR.phi.kalman(x=c(0.8,0.2),y1=y1,y2=y2,t=st,yerr1=yerr1,yerr2=yerr2)
BIAR.phi.kalman<-function (x, y1, y2, t, yerr1, yerr2, zero.mean = "TRUE")
{
  sigmay <- var(cbind(y1,y2))
  if (zero.mean == "FALSE")
  {
    y1 = y1 - mean(y1)
    y2 = y2 - mean(y2)
  }
  n = length(y1)
  Sighat = sigmay %*% matrix(c(1, 0, 0, 1), 2, 2)
  xhat = matrix(0, nrow = 2, ncol = n)
  delta <- diff(t)
  Q = Sighat
  phi.R = x[1]
  phi.I = x[2]
  F = matrix(0, nrow = 2, ncol = 2)
  G = diag(2)
  phi = complex(1, real = phi.R, imaginary = phi.I)
  sum.Lambda = 0
  sum.error = 0
  phi = ifelse(is.na(phi) == TRUE, 1.1, phi)
  y=rbind(y1,y2)
  if (Mod(phi) < 1) {
    Phi = Mod(phi)
    psi <- acos(phi.R/Phi)
    if (phi.I < 0 & Mod(phi) < 1)
      psi = -acos(phi.R/Phi)
    for (i in 1:(n - 1)) {
      R= matrix(c(yerr1[i + 1]^2,0,0,yerr2[i + 1]^2),2,2)
      Lambda = G %*% Sighat %*% t(G) + R
      if (det(Lambda) <= 0 | length(which(is.na(Lambda)))>0) {
        sum.Lambda <- n * 1e+10
        break
      }
      phi2.R <- (Phi^delta[i]) * cos(delta[i] * psi)
      phi2.I <- (Phi^delta[i]) * sin(delta[i] * psi)
      F[1, 1] = phi2.R
      F[1, 2] = -phi2.I
      F[2, 1] = phi2.I
      F[2, 2] = phi2.R
      phi2 <- 1 - Mod(phi^delta[i])^2
      Qt <- phi2 * Q
      sum.Lambda = sum.Lambda + log(det(Lambda))
      Theta=F%*%Sighat%*%t(G)
      sum.error = sum.error + t(y[,i] - G %*% xhat[, i]) %*% solve(Lambda) %*% (y[,i] - G %*% xhat[, i])
      xhat[,i+1]=F%*%xhat[,i]+Theta%*%solve(Lambda)%*%(y[,i]-G%*%xhat[,i])
      if(sum(R)==0)
        Sighat=Qt + F%*%(Sighat-t(Sighat))%*%t(F)
      else
        Sighat=F%*%Sighat%*%t(F)+ Qt - Theta%*%solve(Lambda)%*%t(Theta)
    }
    yhat = G %*% xhat
    out <- ifelse(is.na(sum.Lambda) == TRUE, 1e+10, (sum.Lambda +
                                                       sum.error)/2)
  }
  else out = 1e+10
  return(out)
}
