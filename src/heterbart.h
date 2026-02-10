/*
 */

//#define ARMA_USE_SUPERLU 1
#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h
#include "bart.h"
#include "heterbartfuns.h"
#include "heterbd.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

class heterbart : public bart
{
public:
  heterbart():bart() {}
  heterbart(size_t m):bart(m) {}
  void pr();
  void draw(double *sigma, rn& gen, double nu, double lambda);
  int findindex(Rcpp::IntegerVector vec, int val);
  void updatetrees(Rcpp::List treeprev);
  void updatetrees_module(Rcpp::List treeprev, int rowid, int treeid);
  double draw_sigmaupdate(double *sigma, rn& gen, double nu, double lambda, double &kappa, double &sigma_m, double &mlik, bool isexact);
  void draw_test(double *sigma, rn& gen, double &kappa, double &sigma_m, bool isexact);
  void initializetrees();
  arma::sp_mat convertSparse(Rcpp::S4 mat);
  void makemesh();
  void makemesh_bypoints();
  void makemesh_bypoints_unique();
  Rcpp::List printtrees();
};

#endif
