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

//#define ARMA_USE_SUPERLU
#include "info.h"
#include "makeC.h"
#include "RcppArmadillo.h"
#include "heterbartfuns.h"


//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.

// maybe the first term should be positive?
double heterlh(double b, double M, double tau) {
  double t2 =tau*tau;
  double k = b*t2+1;
  return -.5*log(k)+.5*M*M*t2/k;
}
//--------------------------------------------------
//compute b and M  for left and right give bot and v,c
void hetergetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& bl, double& Ml, size_t& nr,  double& br, double& Mr, double *sigma, double *r)
{
  double *xx;//current x
  bl=0; Ml=0.0; br=0; Mr=0.0; nl=0; nr=0;
  double w;
  // iterate through every observation
  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
      w= 1.0/(sigma[i]*sigma[i]);
      if(xx[v] < xi[v][c]) {
        nl+=1;
        bl+=w;
        //cout << r[i] << " " << di.y[i] << " " << endl;
        //Ml += w*r[i];
        Ml += w*di.y[i];
      } else {
        nr+=1;
        br+=w;
        Mr += w*di.y[i];
        //Mr += w*r[i];
      }
    }
  }
}

//compute b and M  for left and right give bot and v,c
void countbotnodes(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl,  size_t& nr)
{
  double *xx = di.x;//current x
  nl=0; nr=0;
  for(size_t i=0;i<di.n;i++) {
    //xx = di.x + i*di.p;
    if(nx==x.bn(xx+i*di.p,xi)) { //does the bottom node = xx's bottom node
      if(*(xx+i*di.p+v) < xi[v][c]) {
        nl+=1;
      } else {
        nr+=1;
      }
    }
  }
}

//--------------------------------------------------
//compute b and M for left and right bots
void hetergetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, double& bl, double& Ml, double& br, double& Mr, double *sigma)
{

  double *xx;//current x
  bl=0; Ml=0.0; br=0; Mr=0.0;
  double w;

  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    tree::tree_cp bn = x.bn(xx,xi);
    if(bn==l) {
      w = 1.0/(sigma[i]*sigma[i]);
      bl+=w;
      Ml += w*di.y[i];
    }
    if(bn==r) {
      w = 1.0/(sigma[i]*sigma[i]);
      br+=w;
      Mr += w*di.y[i];
    }
  }
}
//--------------------------------------------------
//draw one mu from post
double heterdrawnodemu(double b, double M, double tau, rn& gen)
{
  double muhat = M/b;
  double a = 1.0/(tau*tau);
  return (b*muhat)/(a+b) + gen.normal()/sqrt(a+b);
}
//--------------------------------------------------
// I changed sigma* to sigma so sigma[0] becomes sigma
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<double>& bv, std::vector<double>& Mv,double *sigma)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x

  bnv.clear();
  x.getbots(bnv);

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();
  bv.resize(nb);
  Mv.resize(nb);

  std::map<tree::tree_cp,size_t> bnmap;
  for(bvsz i=0;i!=bnv.size();i++) {bnmap[bnv[i]]=i;bv[i]=0;Mv[i]=0.0;}

  double w;
  for(size_t i=0;i<di.n;i++) {
    w = 1.0/(sigma[i]*sigma[i]);
    xx = di.x + i*di.p;
    tbn = x.bn(xx,xi);
    ni = bnmap[tbn];

    bv[ni] += w;
    Mv[ni] += w*di.y[i];
  }
}
//--------------------------------------------------
//heter version of drmu, need b and M instead of n and sy
void heterdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
  tree::npv bnv;
  std::vector<double> bv;
  std::vector<double> Mv;
  heterallsuff(t,xi,di,bnv,bv,Mv,sigma);
  for(tree::npv::size_type i=0;i!=bnv.size();i++)
    bnv[i]->settheta(heterdrawnodemu(bv[i],Mv[i],pi.tau,gen));
}

Rcpp::DataFrame getCoordsDf(dinfo& di) {
  Rcpp::NumericVector s1 = di.s1;
  Rcpp::NumericVector s2 = di.s2;
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  myList[0] = s1;
  namevec.push_back("cx");
  myList[1] = s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  return(dfout);
}

Rcpp::DataFrame getUniqueCoordsDf(dinfo& di) {
  Rcpp::NumericVector s1 = di.s1;
  Rcpp::NumericVector s2 = di.s2;
  Rcpp::NumericVector s1unique;
  Rcpp::NumericVector s2unique;
  int pos = 0;
  for (int i = 0; i < di.nunique; i++) {
    s1unique.push_back(s1[pos]);
    s2unique.push_back(s2[pos]);
    pos += di.dat_size[i];
  }
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  myList[0] = s1unique;
  namevec.push_back("cx");
  myList[1] = s2unique;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  return(dfout);
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat calcPrecMat(dinfo& di, pinfo& pi, double sigma) {
  //arma::arma_rng::set_seed(1);
  Rcpp::DataFrame dfout = getCoordsDf(di);
  Environment env("package:BARTSIMP");
  Function maternmat = env["matern.mat"];
  Rcpp::NumericMatrix cov = maternmat(Rcpp::_["nu"] = pi.nu,
                                      Rcpp::_["kappa"] = pi.kappa,
                                      Rcpp::_["coords"] = dfout,
                                      Rcpp::_["sigma2e"] = 0,
                                      Rcpp::_["sigma2x"] = pi.sigma2x);
  arma::mat x_ = as<arma::mat>(cov); // Sigma
  arma::mat y_(di.n, di.n, fill::eye);
  y_ = y_* sigma * sigma;
  arma::mat res = (x_ + y_).i();
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat calcCovMat(dinfo& di, pinfo& pi, double sigma) {
  //arma::arma_rng::set_seed(1);
  Rcpp::DataFrame dfout = getCoordsDf(di);
  Environment env("package:BARTSIMP");
  Function maternmat = env["matern.mat"];
  Rcpp::NumericMatrix cov = maternmat(Rcpp::_["nu"] = pi.nu,
                                      Rcpp::_["kappa"] = pi.kappa,
                                      Rcpp::_["coords"] = dfout,
                                      Rcpp::_["sigma2e"] = 0,
                                      Rcpp::_["sigma2x"] = pi.sigma2x);
  arma::mat x_ = as<arma::mat>(cov); // Sigma
  arma::mat y_(di.n, di.n, fill::eye);
  y_ = y_* sigma * sigma;
  arma::mat res = (x_ + y_);
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat calcMat_weighted(dinfo& di, pinfo& pi, double sigma, double sigmax, double kappa) {
  //arma::arma_rng::set_seed(1);
  Rcpp::DataFrame dfout = getUniqueCoordsDf(di);
  Environment env("package:BARTSIMP");
  Function maternmat = env["matern.mat"];
  Rcpp::NumericMatrix cov = maternmat(Rcpp::_["nu"] = pi.nu,
                                      Rcpp::_["kappa"] = kappa,
                                      Rcpp::_["coords"] = dfout,
                                      Rcpp::_["sigma2e"] = 0,
                                      Rcpp::_["sigma2x"] = sigmax*sigmax);
  arma::mat x_ = as<arma::mat>(cov); // Sigma
  arma::mat y_(di.nunique, di.nunique, fill::eye);
  for (int i = 0; i < di.nunique; i++) {
    y_(i,i) = y_(i,i)* sigma * sigma/di.dat_size[i];
  }
  arma::mat res = (x_ + y_);
  return res;
}

arma::mat calcPrecMat_test(dinfo& di, pinfo& pi, double sigma) {
  Rcpp::DataFrame dfout = getCoordsDf(di);
  Environment env("package:BARTSIMP");
  Function maternmat = env["matern.mat"];
  arma::mat y_(di.n, di.n, fill::eye);
  y_ = y_* (sigma * sigma);
  arma::mat res = (y_).i();
  return res;
}

std::vector<double> calcmeans(dinfo& di, double *r) {
  double nl = di.dat_size.size();
  double pos = 0;
  double posend;
  std::vector<double> ybar;
  for (int i = 0; i < nl; i++) {
    posend = pos + di.dat_size[i];
    //cout << pos<< " " << posend << endl;
    double temp = 0;
    for (int j = pos; j < posend; j++) {
      temp += r[j];
    }
    temp = temp/di.dat_size[i];
    ybar.push_back(temp);
    pos += di.dat_size[i];
  }
  return ybar;
}

void printXunique(dinfo& di) {
  double* x = di.xunique;
  std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt", MyFile.out | MyFile.app);
  for (int i = 0; i < di.nunique; i++) {
    for (int j = 0; j < di.p; j++) {
      MyFile << *(x + di.p*i+j) << " ";
    }
    MyFile << endl;
  }
  MyFile.close();
}


void heterdrmu_new(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
                   rn& gen, double *r, double &sigma_m, double &kappa) {
  // convert to std::vector
  //arma::arma_rng::set_seed(1);
  //std::vector<double> yhat_;
  //for (int i = 0; i < di.n; i++) {
  //  yhat_.push_back(r[i]);
  //}
  std::vector<double> yhat_ = calcmeans(di, r);
  arma::vec yhat(yhat_);
  std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt", MyFile.out | MyFile.app);
  // MyFile << "---- Mu update ----" << endl;
  // for (int i = 0; i < yhat.size(); i++) {
  //   MyFile << yhat[i] << " ";
  // }
  // MyFile << endl;
  // create the C matrix
  std::vector<size_t> idv1;
  //printXunique(di);
  getbotsid(t, di.xunique, xi,di.nunique , di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.nunique);
  //printC(df);
  Rcpp::NumericMatrix dfmat = testDFtoNM(df);
  arma::mat Cmat = convertC(dfmat);
  // MyFile << "==== spat var, matern var, kappa ====" << endl;
  // MyFile << sigma << " " << sigma_m  << " " << kappa << endl;
  arma::mat prec = calcMat_weighted(di, pi, sigma, sigma_m, kappa);
  // MyFile << "==== prec ====" << endl;
  // for (int i = 0; i < prec.n_rows; i++) {
  //   for (int j = 0; j < prec.n_cols; j++) {
  //     MyFile << prec(i,j) << " ";
  //   }
  //   MyFile << endl;
  // }
  arma::mat temp1 = Cmat.t()*prec;
  arma::mat temp2(Cmat.n_cols , Cmat.n_cols, fill::eye);
  temp2 = temp2 * pi.tau * pi.tau;
  temp2 = (temp1*Cmat + temp2).i();
  arma::mat L = arma::chol(temp2, "lower");
  //MyFile << "==== L ====" << endl;
  // for (int i = 0; i < L.n_rows; i++) {
  //   for (int j = 0; j < L.n_cols; j++) {
  //     MyFile << L(i,j) << " ";
  //   }
  //   MyFile << endl;
  // }
  arma::vec meanvec = temp2*temp1*yhat;
  //arma::arma_rng::set_seed_random();
  GetRNGstate();
  Rcpp::NumericVector gmu = Rcpp::rnorm(Cmat.n_cols);
  PutRNGstate();
  arma::vec tempz(Cmat.n_cols);
  for (int i = 0; i < tempz.size(); i++) {
    tempz(i) = gmu[i];
  }
  MyFile << "==== tempz ====" << endl;
  for (int i = 0; i < tempz.size(); i++) {
    MyFile << tempz(i) << " ";
  }
  MyFile << endl;
  arma::vec mudraws = meanvec + L*tempz;
  //MyFile << "==== mudraws ====" << endl;
  // for (int i = 0; i < mudraws.size(); i++) {
  //   MyFile << mudraws(i) << " ";
  // }
  // MyFile << endl;
  // get a list of unique bot ids
  std::set<int> s;
  unsigned size = idv1.size();
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv_unique;
  idv_unique.assign(s.begin(),s.end());
  //MyFile << "==== idvunique ====" << endl;
  for (int i = 0; i < idv_unique.size(); i++) {
    //MyFile << idv_unique[i] << " ";
    int btid = idv_unique[i];
    tree::tree_p np = t.getptr(btid);
    np->settheta(mudraws[i]);
  }
  // MyFile << endl;
  // MyFile.close();
  //t.printtreeToFile(xi);
  //arma::mat prec =
  return;
}

// covert to sparse matrix
arma::sp_mat convertSparse_forQ(Rcpp::S4 mat) {
  // obtain dim, i, p. x from S4 object

  IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec j = Rcpp::as<arma::urowvec>(mat.slot("j"));
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

  int nrow = dims[0];
  int ncol = dims[1];
  // use Armadillo sparse matrix constructor
  arma::sp_mat res(i, j, x, nrow, ncol);

  return res;
  // Rcout << "SpMat res:\n" << res << std::endl;
}

void heterdrmu_new_approx(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
                   rn& gen, double *r, double &sigma_m, double &kappa) {
  // convert to std::vector
  //arma::arma_rng::set_seed(1);
  //std::vector<double> yhat_;
  //for (int i = 0; i < di.n; i++) {
  //  yhat_.push_back(r[i]);
  //}
  std::vector<double> yhat_ = calcmeans(di, r);
  arma::vec yhat(yhat_);

  //std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt", MyFile.out | MyFile.app);
  // MyFile << "---- Mu update ----" << endl;
  // for (int i = 0; i < yhat.size(); i++) {
  //   MyFile << yhat[i] << " ";
  // }
  // MyFile << endl;
  // create the C matrix
  std::vector<size_t> idv1;
  //printXunique(di);
  getbotsid(t, di.xunique, xi,di.nunique , di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.nunique);
  printC(df);
  Rcpp::NumericMatrix dfmat = testDFtoNM(df);
  arma::mat Cmat = convertC(dfmat);

  // call INLA from Rcpp
  Environment env("package:INLA");
  // define the pcmatern model
  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];
  Function inlaSpdePrecision= env["inla.spde.precision"];
  // SPDE stuff
  Rcpp::NumericVector priorRange = {sqrt(pi.nu*8.0)/kappa,NA_REAL};
  Rcpp::NumericVector priorSigma = {sigma_m,NA_REAL};
  Rcpp::List spde = inlaSpde2PcMatern(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["prior.range"] = priorRange,
    Rcpp::_["prior.sigma"] = priorSigma
  );
  // Generate Q matrix
  Rcpp::NumericVector theta_vec = {-log(4.0*3.1415*sigma_m*sigma_m*kappa*kappa),
                                        log(kappa)};
  Rcpp::S4 Q = inlaSpdePrecision(
    Rcpp::_["spde"] = spde,
    Rcpp::_["theta"] = theta_vec
  );
  arma::sp_mat Q_sp = Rcpp::as<arma::sp_mat>(Q);
  arma::sp_mat A_sp = di.A_unique;
  arma::sp_mat mat_1 = (1/sigma/sigma*A_sp.t()*A_sp + Q_sp);
  arma::mat mydiag(mat_1.n_rows , mat_1.n_cols, fill::eye);
  arma::mat mat_1_inv = arma::spsolve(mat_1, mydiag, "lapack");
  arma::mat mat_2 = A_sp*mat_1_inv*A_sp.t()/(sigma*sigma*sigma*sigma);
  arma::mat mat_3(A_sp.n_rows, A_sp.n_rows, fill::eye);
  arma::mat mat_temp = mat_3/sigma/sigma - mat_2;
  // calculating V matrix
  arma::mat diag_temp(Cmat.n_cols , Cmat.n_cols, arma::fill::eye);
  // cout << Cmat.n_cols << endl;
  // cout << Cmat.n_rows << endl;
  // cout << pi.tau << endl;
  // cout << mat_temp.n_cols << endl;
  // cout << diag_temp.n_cols << endl;
  arma::mat V_inv = Cmat.t()*mat_temp*Cmat + pi.tau*pi.tau*diag_temp;
  arma::vec m = V_inv*Cmat.t()*mat_temp*yhat;
  arma::mat L = arma::chol(V_inv, "lower");
  arma::vec tempz(Cmat.n_cols, arma::fill::randn);
  arma::vec tempv = solve(L, tempz);
  arma::vec mudraws = m + tempv;


  // MyFile << "==== spat var, matern var, kappa ====" << endl;
  // MyFile << sigma << " " << sigma_m  << " " << kappa << endl;
  // arma::mat prec = calcMat_weighted(di, pi, sigma, sigma_m, kappa);
  // MyFile << "==== prec ====" << endl;
  // for (int i = 0; i < prec.n_rows; i++) {
  //   for (int j = 0; j < prec.n_cols; j++) {
  //     MyFile << prec(i,j) << " ";
  //   }
  //   MyFile << endl;
  // }
  // arma::mat temp1 = Cmat.t()*prec;
  // arma::mat temp2(Cmat.n_cols , Cmat.n_cols, fill::eye);
  // temp2 = temp2 * pi.tau * pi.tau;
  // temp2 = (temp1*Cmat + temp2).i();
  // arma::mat L = arma::chol(temp2, "lower");
  //MyFile << "==== L ====" << endl;
  // for (int i = 0; i < L.n_rows; i++) {
  //   for (int j = 0; j < L.n_cols; j++) {
  //     MyFile << L(i,j) << " ";
  //   }
  //   MyFile << endl;
  // }
  // arma::vec meanvec = temp2*temp1*yhat;
  // arma::arma_rng::set_seed_random();
  // arma::vec tempz(Cmat.n_cols, arma::fill::randn);
  // MyFile << "==== tempz ====" << endl;
  // for (int i = 0; i < tempz.size(); i++) {
  //   MyFile << tempz(i) << " ";
  // }
  // MyFile << endl;
  // arma::vec mudraws = meanvec + L*tempz;
  //MyFile << "==== mudraws ====" << endl;
  // for (int i = 0; i < mudraws.size(); i++) {
  //   MyFile << mudraws(i) << " ";
  // }
  // MyFile << endl;
  // get a list of unique bot ids
  std::set<int> s;
  unsigned size = idv1.size();
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv_unique;
  idv_unique.assign(s.begin(),s.end());
  //MyFile << "==== idvunique ====" << endl;
  for (int i = 0; i < idv_unique.size(); i++) {
    //MyFile << idv_unique[i] << " ";
    int btid = idv_unique[i];
    tree::tree_p np = t.getptr(btid);
    np->settheta(mudraws[i]);
  }
  // MyFile << endl;
  // MyFile.close();
  //t.printtreeToFile(xi);
  //arma::mat prec =
  return;
}


void heterdrmu_new_test(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen, double *r) {
  // convert to std::vector
  //arma::arma_rng::set_seed(1);
  std::vector<double> yhat_;
  for (int i = 0; i < di.n; i++) {
    yhat_.push_back(r[i]);
  }
  arma::vec yhat(yhat_);
  // create the C matrix
  std::vector<size_t> idv1;
  getbotsid(t, di.x, xi, di.n, di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.n);
  Rcpp::NumericMatrix dfmat = testDFtoNM(df);
  arma::mat Cmat = convertC(dfmat);
  arma::mat prec = calcPrecMat_test(di, pi, sigma);
  arma::mat temp1 = Cmat.t()*prec;
  arma::mat temp2(Cmat.n_cols , Cmat.n_cols, fill::eye);
  temp2 = temp2 * pi.tau *pi.tau;
  temp2 = (temp1*Cmat + temp2).i();
  arma::mat L = arma::chol(temp2, "lower");
  arma::vec meanvec = temp2*temp1*yhat;
  //arma::arma_rng::set_seed_random();
  arma::vec tempz(Cmat.n_cols, arma::fill::randn);
  arma::vec mudraws = meanvec + L*tempz;
  // get a list of unique bot ids
  std::set<int> s;
  unsigned size = idv1.size();
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv_unique;
  idv_unique.assign(s.begin(),s.end());
  for (int i = 0; i < idv_unique.size(); i++) {
    int btid = idv_unique[i];
    tree::tree_p np = t.getptr(btid);
    np->settheta(mudraws[i]);
  }
  //arma::mat prec =
  return;
}



