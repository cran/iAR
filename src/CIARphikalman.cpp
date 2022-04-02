#include <RcppArmadillo.h>

#include <string.h>
#include <iostream>
#include <cstdio>
#include <complex>
#include <iomanip>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood of the CIAR Model
//'
//' This function return the negative log likelihood of the CIAR process given specific values of phiR and phiI
//'
//' @param x An array with the parameters of the CIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of CIAR model.
//' @param y Array with the time series observations.
//' @param t Array with the irregular observational times.
//' @param yerr Array with the measurements error standard deviations.
//' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
//' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
//' @param c Nuisance parameter corresponding to the variance of the imaginary part.
//'
//' @return Value of the negative log likelihood evaluated in phiR and phiI.
//' @export
//' @references
//' \insertRef{Elorrieta_2019}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{CIARsample}}
//'
//'
//' @examples
//' n=300
//' set.seed(6714)
//' st<-gentime(n)
//' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
//' y=x$y
//' yerr=rep(0,n)
//' CIARphikalman(x=c(0.8,0),y=y,t=st,yerr=yerr)
// [[Rcpp::export]]
double CIARphikalman(arma::vec x, arma::vec y, arma::vec t, arma::vec yerr, String zeroMean="TRUE", String standarized="TRUE",double c=1.0 ) {
  arma::cx_double phi(x[0], x[1]);

  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return 1e10;
  }

  double sigmaY = 1.0;

  if(standarized == "FALSE") {
    sigmaY = arma::var(y);
  }

  if(zeroMean == "FALSE") {
    y = y - arma::mean(y);
  }
  int n = y.size();
  arma::mat auxSigmaHat = {{1, 0}, {0, c}};
  arma::mat sigmaHat = sigmaY * auxSigmaHat;
  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(t);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G(1, 2, fill::zeros);
  G.col(0) = 1;

  double psi = -acos(phi.real()/phiMod);
  double sumError = 0.0;
  double sumLambda = 0.0;

  for (int i = 0; i < (n-1); ++i) {
    arma::mat Lambda = G * sigmaHat * G.t() + pow(yerr[i+1], 2.0);
    if(arma::as_scalar(Lambda) <= 0 || Lambda.has_nan()) {
      sumLambda = n * 1e10;
      break;
    }

    double deltaPsi = delta[i] * psi;
    double phiModDelta = pow(phiMod, delta[i]);

    double phi2Real = phiModDelta * cos(deltaPsi);
    double phi2Imag = phiModDelta * sin(deltaPsi);

    arma::mat F = {{phi2Real, -phi2Imag}, {phi2Imag, phi2Real}};

    arma::cx_double phi2Inter = pow(phi, delta[i]);
    double phi2Mod = std::pow(std::sqrt(std::pow(phi2Inter.real(), 2.0) + std::pow(phi2Inter.imag(), 2.0)), 2.0);
    double phi2 = 1 - phi2Mod;
    auto Qt = phi2 * Q;

    sumLambda = sumLambda + log(arma::as_scalar(Lambda));

    arma::mat theta = F * sigmaHat * G.t();
    arma::mat aux = y.row(i) - (G * xhat.col(i));
    arma::mat auxLambda = (pow(aux, 2.0))/arma::as_scalar(Lambda);

    sumError = sumError + arma::as_scalar(auxLambda);

    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;
    sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
  }

  if(sumLambda != sumLambda) {
    return 1e10;
  } else {
    return (sumLambda + sumError)/n;
  }
}
