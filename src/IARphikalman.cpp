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
//' @param zeroMean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
//' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.
//' @param yest The estimate of a missing value in the time series. This function recognizes a missing value with a NA. If the time series does not have a missing value, this value does not affect the computation of the likelihood.
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
//' IARphikalman(x=0.8,y=y,yerr=yerr,st=st,yest=0)
// [[Rcpp::export]]
double IARphikalman(arma::vec yest,double x, arma::vec y, arma::vec yerr, arma::vec st, bool zeroMean = true, bool standardized = true) {
  arma::vec y_copy = y;

  //Removing NA elements before join them
  uvec y_NA_indexes = arma::find_nonfinite(y);

  if(y_NA_indexes.size() > 0) {
    y_copy.shed_rows(y_NA_indexes);
  }

  double sigmay = 1.0;

  if(zeroMean == false) {
    y = y - arma::mean(y_copy);
  }


  if(standardized == false) {
    sigmay = arma::var(y_copy);
  }

  int n = y.size();

  arma::vec delta(n, fill::zeros);
  arma::mat xhat(1, n, fill::zeros);
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

      arma::vec yaux = y.row(i);
      arma::mat innov;

      innov = yaux - G * xhat.col(i);

      if(yaux.has_nan() == true) {
        arma::mat temp(1, 1, fill::zeros);
        temp.row(0) = yest;
        innov = temp - G * xhat.col(i);
      }


      sumError = sumError + arma::as_scalar((arma::pow(innov, 2)/Lambda));
      xhat.col(i + 1) = F * xhat.col(i) + Theta * arma::inv(Lambda) * innov;
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
