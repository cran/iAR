#include <RcppArmadillo.h>
#include <string.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;


void uvec_push(arma::uvec & v, unsigned int value) {
  arma::uvec av(1);
  av.at(0) = value;
  v.insert_rows(v.n_rows, av.row(0));
}


//' Full Minus Log Likelihood of the BIAR Model
//'
//' This function return the full negative log likelihood of the BIAR process given specific values of phiR and phiI
//'
//'
//' @param yest An array with the estimate of a missing value in one or both time series of the bivariate process. This function recognizes a missing value with a NA. If the bivariate time series does not have a missing value, this value does not affect the computation of the likelihood.
//' @param phiValues An array with the parameters of the BIAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of BIAR model.
//' @param y1 Array with the observations of the first time series of the BIAR process.
//' @param y2 Array with the observations of the second time series of the BIAR process.
//' @param t Array with the irregular observational times.
//' @param yerr1 Array with the measurements error standard deviations of the first time series of the BIAR process.
//' @param yerr2 Array with the measurements error standard deviations of the second time series of the BIAR process.
//' @param zeroMean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
//'
//' @return Value of the full negative log likelihood evaluated in phiR and phiI.
//' @export
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//'
//' @seealso
//'
//' \code{\link{gentime}}, \code{\link{BIARsample}}
//'
//'
//' @examples
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
//' y=x$y
//' y1=y[1,]
//' y2=y[2,]
//' yerr1=rep(0,n)
//' yerr2=rep(0,n)
//' BIARLL(yest=0,phiValues=c(0.8,0.2),y1=y1,y2=y2,t=st,yerr1=yerr1,yerr2=yerr2)
// [[Rcpp::export]]
double BIARLL(arma::vec yest, arma::vec phiValues,
              arma::vec y1, arma::vec y2, arma::vec t, arma::vec yerr1, arma::vec yerr2, String zeroMean = "TRUE") {
  arma::cx_double phi(phiValues[0], phiValues[1]);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return 1e10;
  }

  arma::vec y1_copy = y1;
  arma::vec y2_copy = y2;

  //Removing NA elements before join them
  uvec y1_NA_indexes = arma::find_nonfinite(y1);
  uvec y2_NA_indexes = arma::find_nonfinite(y2);

  if(y1_NA_indexes.size() > 0 && y2_NA_indexes.size() == 0) {
    y1_copy.shed_rows(y1_NA_indexes);
    y2_copy.shed_rows(y1_NA_indexes);
  }

  if(y1_NA_indexes.size() == 0 && y2_NA_indexes.size() > 0) {
    y1_copy.shed_rows(y2_NA_indexes);
    y2_copy.shed_rows(y2_NA_indexes);
  }

  if(y1_NA_indexes.size() > 0 && y2_NA_indexes.size() > 0) {
    uvec combination = join_cols(y1_NA_indexes, y2_NA_indexes);
    uvec all_joined;
    std::set<int> uniqueSet( combination.begin(), combination.end() );

    for(auto elem: uniqueSet) {
      uvec_push(all_joined, elem);
    }

    y1_copy.shed_rows(all_joined);
    y2_copy.shed_rows(all_joined);
  }

  arma::mat y0 = arma::join_horiz(y1_copy, y2_copy);
  arma::mat sigmaY = arma::cov(y0);

  if(zeroMean == "FALSE") {
    y1 = y1 - arma::mean(y1);
    y2 = y2 - arma::mean(y2);
  }

  int n = y1.size();
  arma::mat sigmaHat = sigmaY * arma::eye(2,2);
  arma::mat xhat(2, n, fill::zeros);
  arma::vec delta = arma::diff(t);

  arma::mat Q = sigmaHat;
  arma::mat F(2, 2, fill::zeros);
  arma::mat G = arma::eye(2,2);

  double psi = (phi.imag() < 0 && phiMod < 1) ? -acos(phi.real()/phiMod) : acos(phi.real()/phiMod);

  arma::mat y = arma::join_vert(y1.t(), y2.t());

  double cte = -0.5 * arma::log_det(sigmaHat).real();
  double sumll = cte;

  for (int i = 0; i < (n-1); ++i) {
    arma::vec yaux = y.col(i + 1);
    arma::mat derr = {{std::pow(yerr1[i+1], 2.0), 0}, {0, std::pow(yerr2[i+1], 2.0)}};

    double deltaPsi = delta[i] * psi;
    double phiModDelta = pow(phiMod, delta[i]);

    double phi2Real = phiModDelta * cos(deltaPsi);
    double phi2Imag = phiModDelta * sin(deltaPsi);

    arma::mat F = {{phi2Real, -phi2Imag}, {phi2Imag, phi2Real}};

    arma::cx_double phi2Inter = pow(phi, delta[i]);
    double phi2Mod = std::pow(std::sqrt(std::pow(phi2Inter.real(), 2.0) + std::pow(phi2Inter.imag(), 2.0)), 2.0);
    auto phi2 = 1 - phi2Mod;

    arma::mat Qt = phi2 * Q;
    arma::mat V(2, 1, fill::zeros);

    V.row(0) = (arma::as_scalar(yerr1[i + 1]) > 0) ? Rcpp::rnorm(1, 0, arma::as_scalar(yerr1[i + 1]))[0] : 0;
    V.row(1) = (arma::as_scalar(yerr2[i + 1]) > 0) ? Rcpp::rnorm(1, 0, arma::as_scalar(yerr2[i + 1]))[0] : 0;

    arma::mat aux(2, 1, fill::zeros);
    arma::mat innov(2, 1, fill::zeros);

    if(yaux.row(0).has_nan() == true && yaux.row(1).has_nan() == true) {
        arma::mat temp(2, 1, fill::zeros);
        temp.row(0) = yest;
        temp.row(1) = yest;

        xhat.col(i + 1) = temp - V;
        innov = temp - G * xhat.col(i + 1);
    } else if(yaux.row(0).has_nan() == true && yaux.row(1).has_nan() == false) {
      xhat.row(0).col(i + 1) = yest - V.row(0);
      xhat.row(1).col(i + 1) = yaux.row(1) - V.row(1);

      aux.row(0) = yest;
      aux.row(1) = yaux.row(1);

      innov = aux - G * xhat.col(i + 1);
    } else if(yaux.row(0).has_nan() == false && yaux.row(1).has_nan() == true) {
      xhat.row(0).col(i + 1) = yaux.row(0) - V.row(0);
      xhat.row(1).col(i + 1) = yest - V.row(1);

      aux.row(0) = yaux.row(0);
      aux.row(1) = yest;
      innov = aux - G * xhat.col(i + 1);
    } else {
      xhat.col(i + 1) = yaux - V;
      innov = yaux - G * xhat.col(i + 1);
    }

    arma::mat innov_tr = xhat.col(i + 1) - F * xhat.col(i);
    arma::mat comp2(1, 1, fill::zeros);

    double valQt, signQt;
    arma::log_det(valQt, signQt, Qt);

    if(arma::det(derr) != 0) {
      double valDerr, signDerr;
      log_det(valDerr, signDerr, derr);

      cte = -0.5 * n * valQt - 0.5 * n * valDerr;
      comp2 = -0.5 * innov_tr.t() * (arma::inv(Qt) * innov_tr) - 0.5 * innov.t() * arma::inv(derr) * innov;
    } else {
      cte = -0.5 * n * valQt;
      comp2 = -0.5 * innov_tr.t() * (arma::inv(Qt) * innov_tr);
    }

    sumll = sumll + cte + arma::as_scalar(comp2);
  }

  return -sumll;
}
