/*
 */

#ifndef GUARD_heterbartfuns_h
#define GUARD_heterbartfuns_h

#include "makeC.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"

//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.
double heterlh(double b, double M, double tau);
//--------------------------------------------------
//compute b and M  for left and right give bot and v,c
void hetergetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& bl, double& Ml, size_t& nr, double& br, double& Mr, double *sigma, double *r);
//--------------------------------------------------
//compute b and M for left and right bots
void hetergetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, double& bl, double& Ml, double& br, double& Mr, double *sigma);
//--------------------------------------------------
//draw one mu from post
double heterdrawnodemu(double b, double M, double tau, rn& gen);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<double>& bv, std::vector<double>& Mv, double *sigma);
//--------------------------------------------------
//heter version of drmu, need b and M instead of n and sy
void heterdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);
//------------------
//compute coordinates and store in a dataframe
//Rcpp::DataFrame getCoordsDf(dinfo& di);
Rcpp::DataFrame getCoordsDf(dinfo& di);
Rcpp::DataFrame getUniqueCoordsDf(dinfo& di);
arma::mat calcCovMat(dinfo& di, pinfo& pi, double sigma);
arma::mat calcPrecMat(dinfo& di, pinfo& pi, double sigma);
arma::mat calcPrecMat_test(dinfo& di, pinfo& pi, double sigma);
arma::mat calcMat_weighted(dinfo& di, pinfo& pi, double sigma, double sigmax, double kappa);
std::vector<double> calcmeans(dinfo& di, double *r);
void printXunique(dinfo& di);
void heterdrmu_new(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
                   rn& gen, double *r, double &sigma_m, double &kappa);
void heterdrmu_new_approx(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
                   rn& gen, double *r, double &sigma_m, double &kappa);
arma::sp_mat convertSparse_forQ(Rcpp::S4 mat);
void heterdrmu_new_test(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen, double *r);
void countbotnodes(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl,  size_t& nr);

#endif
