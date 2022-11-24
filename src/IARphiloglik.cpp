#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood of the IAR Model
//'
//' This function return the negative log likelihood of the IAR Model for a specific value of phi.
//'
//' @param x A given phi coefficient of the IAR model.
//' @param y Array with the time series observations.
//' @param st Array with the irregular observational times.
//' @param delta_input Array with the measurements error standard deviations.
//' @param zeroMean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
//' @param standardized logical; if TRUE, the array y was standardized; if FALSE, y contains the raw data
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
//'
//' set.seed(6714)
//' st<-gentime(n=100)
//' y<-IARsample(phi=0.99,st=st,n=100)
//' y<-y$series
//' IARphiloglik(x=0.8,y=y,st=st,delta_input=c(0))
// [[Rcpp::export]]
double IARphiloglik(double x, arma::vec y, arma::vec st, arma::vec delta_input, bool zeroMean = true, bool standardized = true) {
  double sigma = 1;
  int mu = 0;

  int n = y.size();
  arma::vec delta(n - 1, fill::zeros);

  if(arma::sum(delta_input) != 0) {
    delta = delta_input.rows(1, n-1);
  }

  if(standardized==false) {
    sigma = arma::var(y);
  }

  if(zeroMean == false) {
    mu = arma::mean(y);
  }

  arma::vec d = arma::diff(st);
  arma::vec phi(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    phi[i] = std::pow(x, d[i]);
  }

  arma::vec yhat(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    yhat[i] = mu + phi[i] * (y[i] - mu);
  }

  double cte = (n/2) * std::log(2*datum::pi);

  arma::vec temp0 = sigma * (1 - arma::pow(phi, 2)) + arma::pow(delta, 2);
  arma::vec temp1 = arma::log(temp0);
  arma::vec temp2 = (y.rows(1, n-1) - yhat);
  arma::vec temp3 = arma::pow(temp2, 2);

  return cte + 0.5 * arma::sum(temp1 + temp3/temp0);
}
