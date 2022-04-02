#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Forecast from CIAR model
//'
//' Forecast from models fitted by \code{\link{CIARkalman}}
//'
//' @param phiR Real part of the phi coefficient of CIAR model.
//' @param phiI Imaginary part of the phi coefficient of CIAR model.
//' @param y1 Array with the time series observations.
//' @param st Array with the irregular observational times.
//' @param nAhead The number of steps ahead for forecast is required.
//'
//' @return A list with the following components:
//' \itemize{
//' \item{fitted}{ Fitted values by the CIAR model.}
//' \item{forecast}{ Point Forecasts in the n.ahead times.}
//' \item{Lambda}{ Lambda value estimated by the CIAR model at the last time point.}
//' \item{Sighat}{ Covariance matrix estimated by the CIAR model at the last time point.}
//' }
//' @export
//' @references
//' \insertRef{Elorrieta_2019}{iAR}
//'
//' @seealso
//'
//' \code{\link{CIARsample}}, \code{\link{CIARkalman}}, \code{\link{CIARfit}}
//'
//'
//' @examples
//' #Simulated Data
//' n=100
//' set.seed(6714)
//' st<-gentime(n)
//' x=CIARsample(n=n,phiR=0.9,phiI=0,st=st,c=1)
//' y=x$y
//' y1=y/sd(y)
//' n=length(y1)
//' p=trunc(n*0.99)
//' ytr=y1[1:p]
//' yte=y1[(p+1):n]
//' str=st[1:p]
//' ste=st[(p+1):n]
//' n.ahead=ste-str[p]
//'
//' final<-matrix(0,length(n.ahead),4)
//' ciar=CIARkalman(y=ytr,t=str)
//' forCIAR<-CIARforecast(ciar$phiR,ciar$phiI,ytr,str,nAhead=n.ahead)
// [[Rcpp::export]]
List CIARforecast(double phiR, double phiI, arma::vec y1, arma::vec st, double nAhead) {
  List output;

  arma::cx_double phi(phiR, phiI);
  double phiMod = std::sqrt(std::pow(phi.real(), 2.0) + std::pow(phi.imag(), 2.0));
  if(phiMod >= 1.0) {
    return output;
  }

  arma::vec phiValues(2);
  phiValues[0] = phiR;
  phiValues[1] = phiI;

  Function CIARfit("CIARfit");
  List ciarFitOutput = CIARfit(phiValues, y1, st);

  arma::vec yhat = ciarFitOutput["yhat"];
  arma::mat xhat = ciarFitOutput["xhat"];
  arma::vec Theta = ciarFitOutput["Theta"];
  arma::vec Lambda = ciarFitOutput["Lambda"];
  arma::mat Sighat = ciarFitOutput["Sighat"];
  arma::mat Qt = ciarFitOutput["Qt"];

  int n = yhat.size();
  arma::mat G(1, 2, fill::zeros);
  G.row(0).col(0) = 1;

  double psi = -acos(phi.real()/phiMod);

  double phi2R = pow(phiMod, nAhead) * cos(nAhead * psi);
  double phi2I = pow(phiMod, nAhead) * sin(nAhead * psi);

  double yhat1 = 0;
  arma::mat xhat1(2, 1, fill::zeros);
  double Lambda2 = 0;
  arma::mat sigHat2;

  arma::mat F = {{phi2R, -phi2I}, {phi2I, phi2R}};
  arma::vec temp = y1[n-1] - G * xhat.col(n-1);

  xhat1 = (F * xhat.col(n-1)) + Theta * arma::inv(Lambda) * temp;
  yhat1 = arma::as_scalar(G * xhat1);

  sigHat2 = (F * Sighat * F.t()) + Qt - Theta * arma::inv(Lambda) * Theta.t();
  Lambda2 = arma::as_scalar(G * sigHat2 * G.t());

  output["fitted"] = yhat.t();
  output["forecast"] = yhat1;
  output["Lambda"] = Lambda2;
  output["Sighat"] = sigHat2;

  return output;
}
