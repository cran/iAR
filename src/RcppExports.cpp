// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BIARfit
List BIARfit(arma::vec phiValues, arma::vec y1, arma::vec y2, arma::vec t, arma::vec yerr1, arma::vec yerr2, bool zeroMean);
RcppExport SEXP _iAR_BIARfit(SEXP phiValuesSEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP tSEXP, SEXP yerr1SEXP, SEXP yerr2SEXP, SEXP zeroMeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type phiValues(phiValuesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yerr1(yerr1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yerr2(yerr2SEXP);
    Rcpp::traits::input_parameter< bool >::type zeroMean(zeroMeanSEXP);
    rcpp_result_gen = Rcpp::wrap(BIARfit(phiValues, y1, y2, t, yerr1, yerr2, zeroMean));
    return rcpp_result_gen;
END_RCPP
}
// BIARforecast
List BIARforecast(double phiR, double phiI, arma::vec y1, arma::vec y2, arma::vec t, double tAhead);
RcppExport SEXP _iAR_BIARforecast(SEXP phiRSEXP, SEXP phiISEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP tSEXP, SEXP tAheadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phiR(phiRSEXP);
    Rcpp::traits::input_parameter< double >::type phiI(phiISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type tAhead(tAheadSEXP);
    rcpp_result_gen = Rcpp::wrap(BIARforecast(phiR, phiI, y1, y2, t, tAhead));
    return rcpp_result_gen;
END_RCPP
}
// BIARphikalman
double BIARphikalman(arma::vec yest, arma::vec phiValues, arma::vec y1, arma::vec y2, arma::vec t, arma::vec yerr1, arma::vec yerr2, bool zeroMean);
RcppExport SEXP _iAR_BIARphikalman(SEXP yestSEXP, SEXP phiValuesSEXP, SEXP y1SEXP, SEXP y2SEXP, SEXP tSEXP, SEXP yerr1SEXP, SEXP yerr2SEXP, SEXP zeroMeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yest(yestSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phiValues(phiValuesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yerr1(yerr1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yerr2(yerr2SEXP);
    Rcpp::traits::input_parameter< bool >::type zeroMean(zeroMeanSEXP);
    rcpp_result_gen = Rcpp::wrap(BIARphikalman(yest, phiValues, y1, y2, t, yerr1, yerr2, zeroMean));
    return rcpp_result_gen;
END_RCPP
}
// CIARfit
List CIARfit(arma::vec phiValues, arma::vec y, arma::vec t, bool standardized, double c);
RcppExport SEXP _iAR_CIARfit(SEXP phiValuesSEXP, SEXP ySEXP, SEXP tSEXP, SEXP standardizedSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type phiValues(phiValuesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< bool >::type standardized(standardizedSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(CIARfit(phiValues, y, t, standardized, c));
    return rcpp_result_gen;
END_RCPP
}
// CIARforecast
List CIARforecast(double phiR, double phiI, arma::vec y1, arma::vec st, double tAhead);
RcppExport SEXP _iAR_CIARforecast(SEXP phiRSEXP, SEXP phiISEXP, SEXP y1SEXP, SEXP stSEXP, SEXP tAheadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phiR(phiRSEXP);
    Rcpp::traits::input_parameter< double >::type phiI(phiISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< double >::type tAhead(tAheadSEXP);
    rcpp_result_gen = Rcpp::wrap(CIARforecast(phiR, phiI, y1, st, tAhead));
    return rcpp_result_gen;
END_RCPP
}
// CIARphikalman
double CIARphikalman(arma::vec yest, arma::vec x, arma::vec y, arma::vec t, arma::vec yerr, bool zeroMean, bool standardized, double c);
RcppExport SEXP _iAR_CIARphikalman(SEXP yestSEXP, SEXP xSEXP, SEXP ySEXP, SEXP tSEXP, SEXP yerrSEXP, SEXP zeroMeanSEXP, SEXP standardizedSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yest(yestSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yerr(yerrSEXP);
    Rcpp::traits::input_parameter< bool >::type zeroMean(zeroMeanSEXP);
    Rcpp::traits::input_parameter< bool >::type standardized(standardizedSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(CIARphikalman(yest, x, y, t, yerr, zeroMean, standardized, c));
    return rcpp_result_gen;
END_RCPP
}
// CIARsample
List CIARsample(int n, double phiR, double phiI, arma::vec st, int rho, int c);
RcppExport SEXP _iAR_CIARsample(SEXP nSEXP, SEXP phiRSEXP, SEXP phiISEXP, SEXP stSEXP, SEXP rhoSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type phiR(phiRSEXP);
    Rcpp::traits::input_parameter< double >::type phiI(phiISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< int >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(CIARsample(n, phiR, phiI, st, rho, c));
    return rcpp_result_gen;
END_RCPP
}
// IARgsample
List IARgsample(double phi, arma::vec st, int n, int sigma2, int mu);
RcppExport SEXP _iAR_IARgsample(SEXP phiSEXP, SEXP stSEXP, SEXP nSEXP, SEXP sigma2SEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< int >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(IARgsample(phi, st, n, sigma2, mu));
    return rcpp_result_gen;
END_RCPP
}
// IARphigamma
double IARphigamma(arma::vec yest, arma::vec x_input, arma::vec y, arma::vec st);
RcppExport SEXP _iAR_IARphigamma(SEXP yestSEXP, SEXP x_inputSEXP, SEXP ySEXP, SEXP stSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yest(yestSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_input(x_inputSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    rcpp_result_gen = Rcpp::wrap(IARphigamma(yest, x_input, y, st));
    return rcpp_result_gen;
END_RCPP
}
// IARphikalman
double IARphikalman(arma::vec yest, double x, arma::vec y, arma::vec yerr, arma::vec st, bool zeroMean, bool standardized);
RcppExport SEXP _iAR_IARphikalman(SEXP yestSEXP, SEXP xSEXP, SEXP ySEXP, SEXP yerrSEXP, SEXP stSEXP, SEXP zeroMeanSEXP, SEXP standardizedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yest(yestSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yerr(yerrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< bool >::type zeroMean(zeroMeanSEXP);
    Rcpp::traits::input_parameter< bool >::type standardized(standardizedSEXP);
    rcpp_result_gen = Rcpp::wrap(IARphikalman(yest, x, y, yerr, st, zeroMean, standardized));
    return rcpp_result_gen;
END_RCPP
}
// IARphiloglik
double IARphiloglik(double x, arma::vec y, arma::vec st, arma::vec delta_input, bool zeroMean, bool standardized);
RcppExport SEXP _iAR_IARphiloglik(SEXP xSEXP, SEXP ySEXP, SEXP stSEXP, SEXP delta_inputSEXP, SEXP zeroMeanSEXP, SEXP standardizedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta_input(delta_inputSEXP);
    Rcpp::traits::input_parameter< bool >::type zeroMean(zeroMeanSEXP);
    Rcpp::traits::input_parameter< bool >::type standardized(standardizedSEXP);
    rcpp_result_gen = Rcpp::wrap(IARphiloglik(x, y, st, delta_input, zeroMean, standardized));
    return rcpp_result_gen;
END_RCPP
}
// IARphit
double IARphit(arma::vec yest, arma::vec x, arma::vec y, arma::vec st, double nu);
RcppExport SEXP _iAR_IARphit(SEXP yestSEXP, SEXP xSEXP, SEXP ySEXP, SEXP stSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yest(yestSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(IARphit(yest, x, y, st, nu));
    return rcpp_result_gen;
END_RCPP
}
// IARsample
List IARsample(double phi, arma::vec st, int n);
RcppExport SEXP _iAR_IARsample(SEXP phiSEXP, SEXP stSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type st(stSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(IARsample(phi, st, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iAR_BIARfit", (DL_FUNC) &_iAR_BIARfit, 7},
    {"_iAR_BIARforecast", (DL_FUNC) &_iAR_BIARforecast, 6},
    {"_iAR_BIARphikalman", (DL_FUNC) &_iAR_BIARphikalman, 8},
    {"_iAR_CIARfit", (DL_FUNC) &_iAR_CIARfit, 5},
    {"_iAR_CIARforecast", (DL_FUNC) &_iAR_CIARforecast, 5},
    {"_iAR_CIARphikalman", (DL_FUNC) &_iAR_CIARphikalman, 8},
    {"_iAR_CIARsample", (DL_FUNC) &_iAR_CIARsample, 6},
    {"_iAR_IARgsample", (DL_FUNC) &_iAR_IARgsample, 5},
    {"_iAR_IARphigamma", (DL_FUNC) &_iAR_IARphigamma, 4},
    {"_iAR_IARphikalman", (DL_FUNC) &_iAR_IARphikalman, 7},
    {"_iAR_IARphiloglik", (DL_FUNC) &_iAR_IARphiloglik, 6},
    {"_iAR_IARphit", (DL_FUNC) &_iAR_IARphit, 5},
    {"_iAR_IARsample", (DL_FUNC) &_iAR_IARsample, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_iAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
