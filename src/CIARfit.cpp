#include <RcppArmadillo.h>
#include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Fitted Values of CIAR model
//'
//' Fit a CIAR model to an irregularly observed time series.
//'
//' @param phiValues An array with the parameters of the CIAR model. The elements of the array are, in order, the real and the imaginary part of the phi parameter of the CIAR model.
//' @param y Array with the time series observations.
//' @param t Array with the irregular observational times.
//' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series
//' @param c Nuisance parameter corresponding to the variance of the imaginary part.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{yhat}{ Fitted values of the observable part of CIAR model.}
//' \item{xhat}{ Fitted values of both observable part and imaginary part of CIAR model.}
//' \item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
//' \item{Theta}{ Theta array estimated by the CIAR model at the last time point.}
//' \item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
//' \item{Qt}{ Covariance matrix of the state equation estimated by the CIAR model at the last time point.}
//' }
//' @export
//' @references
//' \insertRef{Elorrieta_2019}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{CIARsample}}, \code{\link{CIARphikalman}},\code{\link{CIARkalman}}
//'
//' @examples
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
//' y=x$y
//' y1=y/sd(y)
//' ciar=CIARkalman(y=y1,t=st)
//' ciar
//' yhat=CIARfit(phiValues=c(ciar$phiR,ciar$phiI),y=y1,t=st)
// [[Rcpp::export]]
List CIARfit(arma::vec phiValues, arma::vec y, arma::vec t, bool standardized=true, double c=1) {
  List output;

  arma::cx_double phi(phiValues[0], phiValues[1]);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  double sigmaY = 1;
  if(standardized == false) {
    sigmaY = arma::var(y);
  }

  int n = y.size();
  arma::mat auxSigmaHat = {{1, 0}, {0, c}};
  arma::mat sigmaHat = sigmaY * auxSigmaHat;

  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(t);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G = arma::eye(1,2);
  G.col(0) = 1;

  double psi = -acos(phi.real()/phiMod);

  arma::mat theta;
  arma::mat Lambda(2, 2, fill::zeros);
  arma::mat Qt;

  for (int i = 0; i < (n-1); ++i) {
    Lambda = G * sigmaHat * G.t();

    if(arma::as_scalar(Lambda) <= 0 || Lambda.has_nan()) {
      break;
    }

    double deltaPsi = delta[i] * psi;
    double phiModDelta = pow(phiMod, delta[i]);

    double phi2Real = phiModDelta * cos(deltaPsi);
    double phi2Imag = phiModDelta * sin(deltaPsi);

    arma::mat F = {{phi2Real, -phi2Imag}, {phi2Imag, phi2Real}};

    arma::cx_double phi2Inter = pow(phi, delta[i]);
    double phi2Mod = std::pow(std::sqrt(std::pow(phi2Inter.real(), 2.0) + std::pow(phi2Inter.imag(), 2.0)), 2.0);
    auto phi2 = 1 - phi2Mod;

    Qt = phi2 * Q;
    theta = F * sigmaHat * G.t();

    arma::mat aux = y[i] - (G * xhat.col(i));
    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;

    sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
  }

  arma::mat yhat = G * xhat;

  output["yhat"] = yhat;
  output["xhat"] = xhat;
  output["Sighat"] = sigmaHat;
  output["Theta"] = theta;
  output["Lambda"] = Lambda;
  output["Qt"] = Qt;

  return output;
}
