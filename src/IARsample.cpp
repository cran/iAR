#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Simulate from an IAR Model
//'
//' Simulates an IAR Time Series Model.
//'
//' @param phi A coefficient of IAR model. A value between 0 and 1
//' @param st Array with observational times.
//' @param n Length of the output time series. A strictly positive integer.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{times}{ Array with observation times.}
//' \item{series}{ Array with simulated IAR data.}
//' }
//' @export
//' @references
//' \insertRef{Eyheramendy_2018}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}
//'
//' @examples
//'
//' set.seed(6714)
//' st<-gentime(n=100)
//' y<-IARsample(phi=0.99,st=st, n=100)
//' y<-y$series
//' plot(st,y,type='l')
// [[Rcpp::export]]
List IARsample(double phi, arma::vec st, int n=100) {
  List output;

  arma::mat Sigma(n, n, fill::zeros);

  for(int i = 0; i < n; ++i) {
    arma::vec d(n-i, fill::zeros);
    d = st[i] - st.rows(i, n-1);
    d = arma::abs(d);

    arma::vec temp(n-i, fill::zeros);
    for(int j = 0; j < (n-i); ++j) {
      temp[j] = std::pow(phi, d[j]);
    }

    Sigma.rows(i, (n-1)).col(i) = temp;
    Sigma.cols(i, (n-1)).row(i) = temp.t();

  }

  Function eigenRfunction("eigen");
  NumericMatrix sigmaRcpp = wrap(Sigma);
  Rcpp::List eigen_results = eigenRfunction(Named("x") = sigmaRcpp, Named("symmetric") = true);
  arma::vec eigenValues = Rcpp::as<arma::vec>(eigen_results[0]);
  arma::mat eigenVectors = Rcpp::as<arma::mat>(eigen_results[1]);

  arma::mat temp = arma::eye(n, n);
  temp.diag() = arma::sqrt(eigenValues);

  arma::mat A = eigenVectors * temp * eigenVectors.t();
  arma::vec e = arma::randn(n);
  arma::vec y = arma::vectorise(A * e);

  output["series"] = y.rows(0, n-1);
  output["times"] = st.rows(0, n-1);

  return output;
}
