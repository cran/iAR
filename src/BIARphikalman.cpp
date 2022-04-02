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

//' Minus Log Likelihood of the BIAR Model
//'
//' This function return the negative log likelihood of the BIAR process given specific values of phiR and phiI
//'
//'
//' @param phiValues An array with the parameters of the BIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of BIAR model.
//' @param y1 Array with the observations of the first time series of the BIAR process.
//' @param y2 Array with the observations of the second time series of the BIAR process.
//' @param t Array with the irregular observational times.
//' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
//' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
//' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
//'
//' @return Value of the negative log likelihood evaluated in phiR and phiI.
//' @export
//'
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//' @seealso
//' \code{\link{gentime}}, \code{\link{BIARsample}}
//'
//'
//' @examples
//' n=300
//' set.seed(6714)
//' st<-gentime(n)
//' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
//' y=x$y
//' y1=y[1,]
//' y2=y[2,]
//' yerr1=rep(0,n)
//' yerr2=rep(0,n)
//' BIARphikalman(phiValues=c(0.8,0.2),y1=y1,y2=y2,t=st,yerr1=yerr1,yerr2=yerr2)
// [[Rcpp::export]]
double BIARphikalman(arma::vec phiValues, arma::vec y1, arma::vec y2, arma::vec t, arma::vec yerr1, arma::vec yerr2, String zeroMean = "TRUE") {
  arma::cx_double phi(phiValues[0], phiValues[1]);

  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return 1e10;
  }

  if(zeroMean == "FALSE") {
    y1 = y1 - arma::mean(y1);
    y2 = y2 - arma::mean(y2);
  }

  int n = y1.size();
  arma::mat sigmaY = arma::cov(arma::join_horiz(y1, y2));
  arma::mat sigmaHat = sigmaY * arma::eye(2,2);
  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(t);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G = arma::eye(2,2);

  double psi = (phi.imag() < 0 && phiMod < 1) ? -acos(phi.real()/phiMod) : acos(phi.real()/phiMod);

  arma::mat y = arma::join_vert(y1.t(), y2.t());

  double sumError = 0.0;
  double sumLambda = 0.0;

  for (int i = 0; i < (n-1); ++i) {
    arma::mat R = {{std::pow(yerr1[i+1], 2.0), 0}, {0, std::pow(yerr2[i+1], 2.0)}};
    arma::mat Lambda(2, 2, fill::zeros);
    Lambda = G * sigmaHat * G.t() + R;

    if(det(Lambda) <= 0 || Lambda.has_nan()) {
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

    auto Qt = (1 - phi2Mod) * Q;

    sumLambda = sumLambda + log(arma::det(Lambda));
    arma::mat theta = F * sigmaHat * G.t();

    arma::mat aux = y.col(i) - (G * xhat.col(i));
    sumError = sumError + as_scalar(aux.t() * arma::inv(Lambda) * aux);

    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;

    if(arma::accu(R) == 0) {
      sigmaHat = Qt + F * (sigmaHat - sigmaHat.t()) * F.t();
    } else {
      sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
    }
  }

  if(sumLambda != sumLambda) {
    return 1e10;
  } else {
    return (sumLambda + sumError)/2;
  }
}
