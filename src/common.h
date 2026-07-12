#ifndef GUARD_common_h
#define GUARD_common_h

#ifdef MATHLIB_STANDALONE
#define NoRcpp
#else
#define RNG_Rcpp
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef NoRcpp

#include <stdio.h> // for printf

using std::cout;

#define PI 3.141592653589793238462643383280

#else // YesRcpp

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#define printf Rprintf
#define cout Rcpp::Rcout

#endif

// log(2*pi)
#define LTPI 1.837877066409345483560659472811
// sqrt(2*pi)
#define RTPI 2.506628274631000502415765284811

#include "rn.h"

#endif
