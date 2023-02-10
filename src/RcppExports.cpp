// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MMCL
List MMCL(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0, const double& tht0, int frailty, int penalty, double tune, int a, int b, int p, double power);
RcppExport SEXP _frailtyMMpen_MMCL(SEXP ySEXP, SEXP XSEXP, SEXP dSEXP, SEXP coef0SEXP, SEXP lambda0SEXP, SEXP tht0SEXP, SEXP frailtySEXP, SEXP penaltySEXP, SEXP tuneSEXP, SEXP aSEXP, SEXP bSEXP, SEXP pSEXP, SEXP powerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type coef0(coef0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< const double& >::type tht0(tht0SEXP);
    Rcpp::traits::input_parameter< int >::type frailty(frailtySEXP);
    Rcpp::traits::input_parameter< int >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< double >::type tune(tuneSEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type power(powerSEXP);
    rcpp_result_gen = Rcpp::wrap(MMCL(y, X, d, coef0, lambda0, tht0, frailty, penalty, tune, a, b, p, power));
    return rcpp_result_gen;
END_RCPP
}
// LogLikCL
double LogLikCL(const NumericVector& y, NumericVector X, const NumericVector& d, const NumericVector& coef0, const NumericVector& lambda0, const double& tht0, int frailty, int a, int b, int p, double power);
RcppExport SEXP _frailtyMMpen_LogLikCL(SEXP ySEXP, SEXP XSEXP, SEXP dSEXP, SEXP coef0SEXP, SEXP lambda0SEXP, SEXP tht0SEXP, SEXP frailtySEXP, SEXP aSEXP, SEXP bSEXP, SEXP pSEXP, SEXP powerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type coef0(coef0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< const double& >::type tht0(tht0SEXP);
    Rcpp::traits::input_parameter< int >::type frailty(frailtySEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type power(powerSEXP);
    rcpp_result_gen = Rcpp::wrap(LogLikCL(y, X, d, coef0, lambda0, tht0, frailty, a, b, p, power));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_frailtyMMpen_MMCL", (DL_FUNC) &_frailtyMMpen_MMCL, 13},
    {"_frailtyMMpen_LogLikCL", (DL_FUNC) &_frailtyMMpen_LogLikCL, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_frailtyMMpen(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
