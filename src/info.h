/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_info_h
#define GUARD_info_h
#include "common.h"
//data
class dinfo {
public:
  dinfo() {p=0;n=0;nunique=0;x=0;y=0;s1=0;s2=0;}
  size_t p;  //number of vars
  size_t n;  //number of observations
  size_t nunique;
  double *x; // jth var of ith obs is *(x + p*i+j)
  double *xunique;
  double *y; // ith y is *(y+i) or y[i]
  Rcpp::NumericVector s1; // x coordinate
  Rcpp::NumericVector s2; // y coordinate
  Rcpp::NumericVector dat_size;
  Rcpp::List mesh; // the mesh object
};
//prior and mcmc
class pinfo
{
public:
  // change rho_0 based on the study range
  // for Kenya level, use rho_0 = 1
  // for unit study region, use 0.1
  pinfo(): pbd(1.0),pb(.5),alpha(.95),mybeta(2.0),tau(1.0),nu(1.0),kappa(sqrt(8)/2.5),
  sigma2e(0.01), sigma2x(0.5), sigmax(sqrt(0.5)), rho_0(2.4), sigma_m0(1.0), alpha_1(0.5), alpha_2(0.5) {}
  //mcmc info
  double pbd; //prob of birth/death
  double pb;  //prob of birth
  //prior info
  double alpha;
  double mybeta;
  double tau;
  //prior about sigma
  //prior for the variance parameters;
  double nu;
  double kappa;
  double sigma2e;
  double sigma2x;
  double sigmax;
  double sigma_m0;
  double rho_0;
  double alpha_1;
  double alpha_2;
  void pr() {
    cout << "pbd,pb: " << pbd << ", " << pb << std::endl;
    cout << "alpha,beta,tau: " << alpha <<
      ", " << mybeta << ", " << tau << std::endl;
  }
};

#endif
