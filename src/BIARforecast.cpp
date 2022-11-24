#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Forecast from BIAR model
//'
//' Forecast from models fitted by \code{\link{BIARkalman}}
//'
//' @param phiR Autocorrelation coefficient of BIAR model.
//' @param phiI Cross-correlation coefficient of BIAR model.
//' @param y1 Array with the observations of the first time series of the BIAR process.
//' @param y2 Array with the observations of the second time series of the BIAR process.
//' @param t Array with the observational times.
//' @param tAhead The time ahead for which the forecast is required.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{fitted}{ Fitted values by the BIAR model.}
//' \item{forecast}{ Point forecast in the time ahead required.}
//' \item{Lambda}{ Lambda value estimated by the BIAR model at the last time point.}
//' \item{Sighat}{ Covariance matrix estimated by the BIAR model at the last time point.}
//' }
//' @export
//' @references
//' \insertRef{Elorrieta_2021}{iAR}
//'
//' @seealso
//'
//' \code{\link{BIARsample}}, \code{\link{BIARkalman}}, \code{\link{BIARfit}}
//'
//'
//' @examples
//' #Simulated Data
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' x=BIARsample(n=n,phiR=0.9,phiI=0.3,st=st)
//' biar=iAR::BIARkalman(y1=x$y[1,],y2=x$y[2,],t=st)
//' forBIAR<-BIARforecast(phiR=biar$phiR,phiI=biar$phiI,y1=x$y[1,],y2=x$y[2,],t=st,tAhead=c(1.3))
// [[Rcpp::export]]
List BIARforecast(double phiR, double phiI, arma::vec y1, arma::vec y2, arma::vec t, double tAhead) {
  List output;

  arma::cx_double phi(phiR, phiI);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  arma::vec phiValues(2);
  phiValues[0] = phiR;
  phiValues[1] = phiI;
  int n = y1.size();
  arma::mat yerr1(1, n, fill::zeros);
  arma::mat yerr2(1, n, fill::zeros);

  Function BIARfit("BIARfit");
  List biarFitOutput = BIARfit(phiValues, y1, y2, t, yerr1, yerr2);

  arma::mat yhat = biarFitOutput["fitted"];
  arma::mat xhat = biarFitOutput["fitted.state"];
  arma::mat theta = biarFitOutput["Theta"];
  arma::mat Lambda = biarFitOutput["Lambda"];
  arma::mat Sighat = biarFitOutput["Sighat"];
  arma::mat Qt = biarFitOutput["Qt"];


  arma::mat G = arma::eye(2,2);
  arma::mat y = arma::join_vert(y1.t(), y2.t());

  arma::mat xhat1(2, 1, fill::zeros);
  arma::mat yhat1(2, 1, fill::zeros);
  arma::mat sigHat2;
  arma::mat Lambda2;

  double psi = (phi.imag() < 0 && phiMod < 1) ? -acos(phi.real()/phiMod) : acos(phi.real()/phiMod);

  double phi2R = pow(phiMod, tAhead) * cos(tAhead * psi);
  double phi2I = pow(phiMod, tAhead) * sin(tAhead * psi);

  arma::mat F = {{phi2R, -phi2I}, {phi2I, phi2R}};
  arma::mat temp = y.col(n-1) - (G * xhat.col(n-1));
  xhat1 = F * xhat.col(n-1) + theta * arma::inv(Lambda) * temp;
  yhat1 = G * xhat1;

  sigHat2 = F * Sighat * F.t() + Qt - theta * arma::inv(Lambda) * theta.t();  
  Lambda2 = G * sigHat2 * G.t();

  output["fitted"] = yhat.t();
  output["forecast"] = yhat1;
  output["Lambda"] = Lambda2;
  output["Sighat"] = sigHat2;

  return output;
}
