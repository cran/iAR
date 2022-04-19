#include <RcppArmadillo.h>
#include <string.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood of the IAR Model estimated via Kalman Recursions
//'
//' This function return the negative log likelihood of the IAR process given a specific value of phi.
//'
//' @param x A given phi coefficient of the IAR model.
//' @param y Array with the time series observations.
//' @param yerr Array with the measurements error standard deviations.
//' @param st Array with the irregular observational times.
//' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
//' @param standarized logical; if true, the array y is standarized; if false, y contains the raw time series.
//'
//' @return Value of the negative log likelihood evaluated in phi.
//' @export
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{IARsample}}
//'
//' @examples
//' set.seed(6714)
//' st<-gentime(n=100)
//' y<-IARsample(phi=0.99,st=st,n=100)
//' y<-y$series
//' yerr=rep(0,100)
//' IARphikalman(x=0.8,y=y,yerr=yerr,st=st)
// [[Rcpp::export]]
double IARphikalman(double x, arma::vec y, arma::vec yerr, arma::vec st, String zeroMean = "FALSE", String standarized = "TRUE") {
  int n = y.size();

  arma::vec delta(n, fill::zeros);
  arma::mat xhat(1, n, fill::zeros);

  double sigmay = 1;
  if(standarized == "FALSE") {
    sigmay = arma::var(y);
  }

  if(zeroMean == "TRUE") {
    y = y - arma::mean(y);
  }

  arma::vec Sighat(1, 1, fill::zeros);

  Sighat[0] = sigmay;
  delta = arma::diff(st);

  arma::vec Q(1, 1, fill::zeros);
  Q[0] = Sighat[0];

  arma::vec F(1,1, fill::zeros);
  arma::vec G(1,1, fill::ones);

  double sumLambda = 0;
  double sumError = 0;

  //Check if x is NaN
  if(x != x) {
    x = 1.1;
  }

  if(std::abs(x) < 1.0) {
    for (int i = 0; i < (n-1); ++i) {
      arma::vec Lambda(1,1, fill::ones);
      arma::vec Theta(1,1, fill::ones);

      Lambda = G * Sighat * G.t() + arma::pow(yerr.row(i+1), 2);

      if(Lambda[0] <= 0 || Lambda.has_nan()) {
        sumLambda = n * 1e10;
        break;
      }

      double phi2 = std::pow(x, arma::as_scalar(delta.row(i)));
      F[0] = phi2;

      phi2 = 1 - (std::pow(x, (2 * arma::as_scalar(delta.row(i)))));
      auto Qt = phi2 * Q;

      sumLambda = sumLambda + arma::as_scalar(arma::log(Lambda));
      Theta = F * Sighat * G.t();

      arma::vec aux = y[i] - (G % xhat.col(i));

      sumError = sumError + arma::as_scalar((arma::pow(aux, 2)/Lambda));
      xhat.col(i + 1) = F * xhat.col(i) + Theta * arma::inv(Lambda) * aux;
      Sighat = F * Sighat * F.t() + Qt - Theta * arma::inv(Lambda) * Theta.t();
    }

    if(sumLambda != sumLambda) {
      return 1e10;
    } else {
      return ((sumLambda + sumError)/n);
    }

  } else {
    return 1e10;
  }
}
