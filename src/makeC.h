#ifndef GUARD_makeC_h
#define GUARD_makeC_h
//#define ARMA_USE_SUPERLU 1

#include "common.h"
//#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// // [[Rcpp::depends(RcppArmadillo)]]

SEXP makeC(std::vector<size_t> idv, size_t n);

void printC(Rcpp::DataFrame df);

NumericMatrix testDFtoNM(DataFrame x);

arma::mat convertC(SEXP mat);

void callPrint(RObject x);

void useOperatorOnVector(NumericVector x);

void useOperatorOnMatrix(NumericMatrix x);

#endif


