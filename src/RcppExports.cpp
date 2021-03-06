// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// GaussianKernel
NumericMatrix GaussianKernel(NumericVector x1, NumericVector x2, List kernelPar);
RcppExport SEXP TempShift_GaussianKernel(SEXP x1SEXP, SEXP x2SEXP, SEXP kernelParSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< List >::type kernelPar(kernelParSEXP);
    __result = Rcpp::wrap(GaussianKernel(x1, x2, kernelPar));
    return __result;
END_RCPP
}
