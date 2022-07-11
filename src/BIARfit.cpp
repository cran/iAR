#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Fitted Values of BIAR model
//'
//' Fit a BIAR model to a bivariate irregularly observed time series.
//'
//' @param phiValues An array with the parameters of the BIAR model. The elements of the array are, in order, the autocorrelation and the cross correlation parameter of the BIAR model.
//' @param y1 Array with the observations of the first time series of the BIAR process.
//' @param y2 Array with the observations of the second time series of the BIAR process.
//' @param t Array with the irregular observational times.
//' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
//' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
//' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{rho}{ Estimated value of the contemporary correlation coefficient.}
//' \item{innov.var}{ Estimated value of the innovation variance.}
//' \item{fitted}{ Fitted values of the BIAR model.}
//' \item{fitted.state}{ Fitted state values of the BIAR model.}
//' \item{Lambda}{ Lambda value estimated by the BIAR model at the last time point.}
//' \item{Theta}{ Theta array estimated by the BIAR model at the last time point.}
//' \item{Sighat}{ Covariance matrix estimated by the BIAR model at the last time point.}
//' \item{Qt}{ Covariance matrix of the state equation estimated by the BIAR model at the last time point.}
//' }
//' @export
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//' @seealso
//' \code{\link{gentime}}, \code{\link{BIARsample}}, \code{\link{BIARphikalman}}, \code{\link{BIARkalman}}
//'
//' @examples
//' \donttest{
//' n=80
//' set.seed(6714)
//' st<-gentime(n)
//' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st,rho=0.9)
//' y=x$y
//' y1=y/apply(y,1,sd)
//' yerr1=rep(0,n)
//' yerr2=rep(0,n)
//' biar=BIARkalman(y1=y1[1,],y2=y1[2,],t=st,delta1 = yerr1,delta2=yerr2)
//' biar
//' predbiar=BIARfit(phiValues=c(biar$phiR,biar$phiI),y1=y1[1,],y2=y1[2,],t=st,yerr1
//'  = rep(0,length(y[1,])),yerr2=rep(0,length(y[1,])))
//' rho=predbiar$rho
//' print(rho)
//' yhat=predbiar$fitted
//' }
// [[Rcpp::export]]
List BIARfit(arma::vec phiValues, arma::vec y1, arma::vec y2, arma::vec t, arma::vec yerr1, arma::vec yerr2, String zeroMean = "TRUE") {
  List output;
  arma::cx_double phi(phiValues[0], phiValues[1]);

  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
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

  arma::mat theta;
  arma::mat Lambda(2, 2, fill::zeros);
  arma::mat Qt;

  for (int i = 0; i < (n-1); ++i) {
    arma::mat R = {{std::pow(yerr1[i+1], 2.0), 0}, {0, std::pow(yerr2[i+1], 2.0)}};
    Lambda = G * sigmaHat * G.t() + R;

    if(det(Lambda) <= 0 || Lambda.has_nan()) {
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

    arma::mat M = phi2 * arma::eye(2,2);
    Qt = M * (Q - R);
    theta = F * sigmaHat * G.t();
    arma::mat aux = y.col(i) - (G * xhat.col(i));

    xhat.col(i + 1) = F * xhat.col(i) + theta * arma::inv(Lambda) * aux;

    if(arma::accu(R) == 0) {
      sigmaHat = Qt + F * (sigmaHat - sigmaHat.t()) * F.t();
    } else {
      sigmaHat = F * sigmaHat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();
    }
  }

  arma::mat yhat = G * xhat;

  arma::mat finalMat = y - G * xhat;
  arma::mat finalCor = arma::cor(finalMat.t());
  arma::mat finalCov = arma::cov(finalMat.t());

  output["rho"] = finalCor.at(0,1);
  output["innov.var"] = finalCov;
  output["fitted"] = yhat;
  output["fitted.state"] = xhat;
  output["Sighat"] = sigmaHat;
  output["Theta"] = theta;
  output["Lambda"] = Lambda;
  output["Qt"] = Qt;

  return output;
}
