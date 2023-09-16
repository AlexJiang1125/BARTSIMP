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

#ifndef GUARD_heterbd_h
#define GUARD_heterbd_h

#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "heterbartfuns.h"

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma,
             std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,double *r);

bool heterbd_test(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double *sigma,
                  std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,
                  double &sigma_m, double &kappa);
bool heterbd_new_test2(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double sigma,
                       std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);
void printX(dinfo& di);
double hetergetmargprob_testonly(
                tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi
);
double hetergetmargprob(tree& x, xinfo& xi, dinfo& di,  double* r, double sigma, pinfo& pi);
double hetergetmargprob_test(
                tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi
);
double hetergetmargprob_all(
                tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi,
                double &sigma_m, double &kappa
);
bool heterbd_new(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double sigma,
                 std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);
bool heterbd_new_test(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double sigma,
                 std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen);
double heterbd_drawsigma_margprob(double* r, dinfo& di, pinfo& pi, double sigma, double &sigma_m, double &kappa);
double heterbd_drawsigma(dinfo& di, double* r, double sigma, double &kappa,double &sigma_m, pinfo& pi, double lambda, double nu);
double heterbd_drawsigma_margprob_test(double* r, dinfo& di, pinfo& pi, double sigma);
double mlikcalcweighted (double* r, dinfo& di, pinfo& pi, double sigma, double &sigma_m, double &kappa);
void heterbd_drawspathyperpars(dinfo& di, double* r, double sigma, double &kappa, double &sigma_m, pinfo& pi, double &mlik);
void heterbd_trihyperpars(dinfo& di, double* r, double& sigma, double &kappa, double &sigma_m, pinfo& pi, double &mlik, double lambda, double nu);
double dmvnorm_distance( arma::vec y, arma::mat Sigma );
double dmvnorm_rcpp( arma::vec y, arma::mat Sigma );

#endif
