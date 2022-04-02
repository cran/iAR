#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;

//' Simulate from an IAR-Gamma Model
//'
//' Simulates an IAR-Gamma Time Series Model.
//'
//' @param phi A coefficient of IAR-Gamma model. A value between 0 and 1.
//' @param st Array with observational times.
//' @param n Length of the output time series. A strictly positive integer.
//' @param sigma2 Scale parameter of the IAR-Gamma process. A positive value.
//' @param mu Level parameter of the IAR-Gamma process. A positive value.
//'
//' @return  A list with the following components:
//' \itemize{
//' \item{y}{ Array with simulated IAR-Gamma process.}
//' \item{st}{ Array with observation times.}
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
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' y<-IARgsample(phi=0.9,st=st,n=n,sigma2=1,mu=1)
//' plot(st,y$y,type='l')
//' hist(y$y,breaks=20)
// [[Rcpp::export]]
List IARgsample(double phi, arma::vec st, int n=100, int sigma2 = 1, int mu = 1) {
  List output;

  arma::vec y(n, fill::zeros);
  arma::vec delta = arma::diff(st);

  y.row(0) = arma::randg(1, distr_param(1,1));

  arma::vec shape(n, fill::zeros);
  arma::vec scale(n, fill::zeros);
  arma::vec yhat(n, fill::zeros);

  for(int i = 1; i < n; ++i) {
    double phid = std::pow(phi, arma::as_scalar(delta.row(i-1)));
    yhat.row(i) = mu + phid * y.row(i-1);
    double gL = sigma2 * (1 - std::pow(phid, 2));
    shape.row(i) = std::pow(arma::as_scalar(yhat.row(i)), 2)/gL;
    scale.row(i) = (gL/yhat.row(i));

    double a = arma::as_scalar(shape.row(i));
    double b = arma::as_scalar(scale.row(i));

    y.row(i) = arma::randg(1, distr_param(a,b));
  }

  output["y"] = y.rows(0, n-1);
  output["st"] = st.rows(0, n-1);

  return output;
}
