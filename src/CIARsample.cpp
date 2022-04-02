#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Simulate from a CIAR Model
//'
//' Simulates a CIAR Time Series Model
//'
//' @param n Length of the output time series. A strictly positive integer.
//' @param st Array with observational times.
//' @param phiR Real part of the coefficient of CIAR model. A value between -1 and 1.
//' @param phiI Imaginary part of the coefficient of CIAR model. A value between -1 and 1.
//' @param rho Correlation between the real and the imaginary part of the process. A value between -1 and 1.
//' @param c Nuisance parameter corresponding to the variance of the imaginary part.
//'
//' @details The chosen phiR and phiI values must satisfy the condition $|phiR + i phiI| < 1$.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{y}{Array with the simulated real part of the CIAR process.}
//' \item{t}{ Array with observation times.}
//' \item{Sigma}{ Covariance matrix of the process.}
//' }
//' @export
//' @references
//' \insertRef{Elorrieta_2019}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}
//'
//' @examples
//' n=300
//' set.seed(6714)
//' st<-gentime(n)
//' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
//' plot(st,x$y,type='l')
//' x=CIARsample(n=n,phiR=-0.9,phiI=0,st=st,c=1)
//' plot(st,x$y,type='l')
// [[Rcpp::export]]
List CIARsample(int n, double phiR, double phiI, arma::vec st, int rho = 0, int c = 1) {
  List output;

  arma::vec delta = arma::diff(st);
  arma::mat x(2, n, fill::zeros);
  arma::mat F(2, 2, fill::zeros);
  arma::mat A(2, 2, fill::zeros);

  arma::cx_double phi(phiR, phiI);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  double psi = -acos(phi.real()/phiMod);

  arma::vec eR = arma::randn(n);
  arma::vec eI = arma::randn(n);

  arma::mat stateError = arma::join_cols(eR.t(), eI.t());

  arma::mat Sigma = {{1, (rho * sqrt(1) * sqrt(c))}, {(rho * sqrt(1) * sqrt(c)), (double) c}};

  arma::mat BU, BV;
  arma::vec Bs;
  arma::svd(BU, Bs, BV, Sigma, "standard");
  A.diag() = Bs;

  arma::mat sigmaRoot = BU * A * BU;

  stateError = stateError.t() * sigmaRoot;
  stateError = stateError.t();

  arma::mat G(1, 2, fill::zeros);
  G.row(0).col(0) = 1;

  arma::vec y(n);
  x.col(0) = stateError.col(0);

  for (int i = 0; i < (n-1); ++i) {
    double deltaPsi = delta[i] * psi;
    double phiModDelta = pow(phiMod, delta[i]);

    double phi2Real = phiModDelta * cos(deltaPsi);
    double phi2Imag = phiModDelta * sin(deltaPsi);
    arma::cx_double phi2Inter = pow(phi, delta[i]);
    double phi2Mod = std::pow(std::sqrt(std::pow(phi2Inter.real(), 2.0) + std::pow(phi2Inter.imag(), 2.0)), 2.0);

    arma::mat F = {{phi2Real, -phi2Imag}, {phi2Imag, phi2Real}};
    double phi2 = 1 - phi2Mod;
    x.col(i + 1) = (F * x.col(i)) + (sqrt(phi2) * stateError.col(i));
    y[i] = arma::as_scalar(G * x.col(i));
  }

  y[n-1] = arma::as_scalar(G * x.col(n-1));

  output["t"] = st;
  output["y"] = y;
  output["sigma"] = Sigma;

  return output;
}
