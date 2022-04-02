#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Minus Log Likelihood IAR-T Model
//'
//' This function return the negative log likelihood of the IAR-T given specific values of phi and sigma.
//'
//' @param x An array with the parameters of the IAR-T model. The first element of the array corresponding to the phi parameter and the second element to the scale parameter sigma
//' @param y Array with the time series observations
//' @param st Array with the irregular observational times
//' @param nu degrees of freedom
//'
//' @return Value of the negative log likelihood evaluated in phi,sigma and nu.
//' @export
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{IARgsample}}
//'
//' @examples
//' n=300
//' set.seed(6714)
//' st<-gentime(n) #Unequally spaced times
//' y<-IARtsample(n,0.9,st,sigma2=1,nu=3)
//' IARphit(x=c(0.9,1),y=y$y,st=st)
// [[Rcpp::export]]
double IARphit(arma::vec x, arma::vec y, arma::vec st, double nu = 3) {
  double sigma = arma::as_scalar(x.row(1));
  int n = y.size();
  arma::vec d = arma::diff(st);

  arma::vec xd(n-1, fill::zeros);
  for(int i = 0; i < n-1; ++i) {
    xd[i] = std::pow(arma::as_scalar(x.row(0)), d[i]);
  }

  arma::vec yhat(n-1, fill::zeros);
  yhat = xd % y.rows(0, n-2);

  arma::vec gL(n-1, fill::zeros);
  gL = sigma * (1 - (arma::pow(xd, 2))) * ((nu-2)/nu);

  double temp1 = std::tgamma((nu+1)/2);
  double temp2 = std::tgamma(nu/2)*(std::sqrt(nu*datum::pi));
  double cte = (n-1) * std::log(temp1/temp2);

  arma::vec stand(n-1, fill::zeros);
  stand = (y.rows(1, n-1) - yhat)/(arma::sqrt(gL));
  stand = arma::pow(stand, 2);

  double s1 = arma::sum(0.5 * arma::log(gL));
  double s2 = arma::sum(arma::log(1 + (1/nu)*stand));

  double out = cte - s1 - ((nu+1)/2) * s2 - 0.5 * (std::log(2*datum::pi) + std::pow(arma::as_scalar(y.row(0)), 2));

  return -out;
}
