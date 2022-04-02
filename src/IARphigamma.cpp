#include <RcppArmadillo.h>
#include <string.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood IAR-Gamma Model
//'
//' This function return the negative log likelihood of the IAR-Gamma given specific values of phi, mu and sigma.
//'
//' @param x_input An array with the parameters of the IAR-Gamma model. The first element of the array corresponding to the phi parameter, the second to the level parameter mu, and the last one to the scale parameter sigma.
//' @param y Array with the time series observations.
//' @param st Array with the irregular observational times.
//'
//' @return Value of the negative log likelihood evaluated in phi, mu and sigma.
//' @export
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{IARgsample}}
//'
//' @examples
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
//' IARphigamma(x_input=c(0.9,1,1),y=y$y,st=st)
// [[Rcpp::export]]
double IARphigamma(arma::vec x_input, arma::vec y, arma::vec st) {
  double x = x_input[0];
  double mu = x_input[1];
  double sigma = x_input[2];

  int n = y.size();
  arma::vec d = arma::diff(st);

  arma::vec xd(st.size(), fill::zeros);
  for(int i = 0; i < st.size(); ++i) {
    xd[i] = std::pow(x, d[i]);
  }

  arma::vec yhat(n-1, fill::zeros);
  for(int i = 0; i < (n-1); ++i) {
    yhat[i] = mu + xd[i] * y[i];
  }

  arma::vec gL(n-1, fill::zeros);
  for(int i = 0; i < (n-1); ++i) {
    gL[i] = sigma * (1 - std::pow(xd[i],2));
  }

  arma::vec beta(n-1, fill::zeros);
  beta=gL/yhat;

  arma::vec alpha(n-1, fill::zeros);
  for(int i = 0; i < (n-1); ++i) {
    alpha[i] = std::pow(yhat[i],2)/gL[i];
  }

  arma::vec temp0 = (-alpha)%arma::log(beta);
  arma::vec temp1 = (alpha - 1) % arma::log(y.rows(1,n-1));

  double out = arma::sum(temp0 - arma::lgamma(alpha) - (y.rows(1,n-1)/beta) + temp1) - arma::as_scalar(y.row(0));
  return -out;
}
