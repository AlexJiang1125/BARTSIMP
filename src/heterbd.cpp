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


#include "heterbd.h"
#include "makeC.h"
#include <iostream>
#include <fstream>

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma,
             std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,
             double *r)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

  if(gen.uniform() < PBx) { //do birth or death

    //--------------------------------------------------
    //draw proposal
    tree::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    //tree xold = x;
    //x.birthp(nx,v,c,0,0);
    //tree xnew = x;
    //x = xold;
    //xnew.printtree(xi);
    //x.printtree(xi);


    //--------------------------------------------------
    //compute sufficient statistics
    size_t nr,nl; //counts in proposed bots
    double bl,br; //sums of weights
    double Ml, Mr; //weighted sum of y in proposed bots
    hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma,r);
    //--------------------------------------------------
    //compute alpha
    //cout << "sigma: " << sigma[0] << " tau: " << pi.tau << endl;
    double alpha=0.0, lalpha=0.0;
    double lhl, lhr, lht;
    if((nl>=5) && (nr>=5)) { //cludge?
      lhl = heterlh(bl,Ml,pi.tau);
      lhr = heterlh(br,Mr,pi.tau);
      lht = heterlh(bl+br,Ml+Mr,pi.tau);
      //cout << "BART:" << nl << " " << nr << endl;
      //double lh1 = hetergetmargprob_test(xnew, xi, di, r, sigma[0], pi);
      //double lh2 = hetergetmargprob_test(x, xi, di, r, sigma[0], pi);
      alpha=1.0;
      lalpha = log(pr) + (lhl+lhr-lht);
      lalpha = std::min(0.0,lalpha);
      //cout << "loglh:" << (lhl+lhr-lht);
      //cout << "loglh2:" << (lh1 - lh2);
    }

    //--------------------------------------------------
    //try metrop
    double mul,mur; //means for new bottom nodes, left and right
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      mul = heterdrawnodemu(bl,Ml,pi.tau,gen);
      mur = heterdrawnodemu(br,Mr,pi.tau,gen);
      x.birthp(nx,v,c,mul,mur);
      //x.birthp(nx,v,c,0,0);
      //cout << nl << " " << nr << endl;
      //cout << "new birth" << endl;
      //x.printtree(xi);
      nv[v]++;
      return true;
    } else {
      //cout << "no birth" << endl;
      //x.printtree(xi);
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal

    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    //tree xold = x;
    //x.deathp(nx,0);
    //tree xnew = x;
    //x = xold;
    //cout << "new" << endl;
    //xnew.printtree(xi);
    //cout << "old" << endl;
    //x.printtree(xi);


    //--------------------------------------------------
    //compute sufficient statistics
    double br,bl; //sums of weights
    double Ml, Mr; //weighted sums of y
    hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma);

    //--------------------------------------------------
    //compute alpha
    double lhl, lhr, lht;
    lhl = heterlh(bl,Ml,pi.tau);
    lhr = heterlh(br,Mr,pi.tau);
    lht = heterlh(bl+br,Ml+Mr,pi.tau);
    //cout << "BART:" << bl*sigma[0]*sigma[0] << " " << br*sigma[0]*sigma[0] << endl;
    //double lh1 = hetergetmargprob_test(xnew, xi, di, r, sigma[0], pi);
    //double lh2 = hetergetmargprob_test(x, xi, di, r, sigma[0], pi);
    double lalpha = log(pr) + (lht - lhl - lhr);
    lalpha = std::min(0.0,lalpha);
    //cout << "loglh:" << -(lhl+lhr-lht);
    //cout << "loglh2:" << (lh1 - lh2);

    //--------------------------------------------------
    //try metrop
    //double a,b,s2,yb;
    double mu;
    if(log(gen.uniform()) < lalpha) {
      mu = heterdrawnodemu(bl+br,Ml+Mr,pi.tau,gen);
      nv[nx->getv()]--;
      x.deathp(nx,mu);
      //x.deathp(nx,0);
      //cout << "new death" << endl;
      //x.printtree(xi);
      return true;
    } else {
      //cout << "no death" << endl;
      //x.printtree(xi);
      return false;
    }
  }
}



// for test only, gaussian case
bool heterbd_test(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double *sigma,
                  std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,
                  double &sigma_m, double &kappa)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

  if(gen.uniform() < PBx) { //do birth or death

    //--------------------------------------------------
    //draw proposal
    tree::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    tree xold = x;
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    size_t nr,nl; //counts in proposed bots
    //compute sufficient statistics
    //size_t nr,nl; //counts in proposed bots
    countbotnodes(x, nx, v, c, xi, di, nl, nr);
    x.birthp(nx,v,c,0,0);
    //bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    //--------------------------------------------------
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0;
    //cout << "nl: " << nl << " " << nr << endl;
    if((nl>=5) && (nr>=5)) { //cludge?
      double lh1 = hetergetmargprob_all(x, xi, di, r, sigma[0], pi, sigma_m, kappa);
      double lh2 = hetergetmargprob_all(xold, xi, di, r, sigma[0], pi, sigma_m, kappa);
      alpha=1.0;
      lalpha = log(pr) + (lh1 - lh2);
      lalpha = std::min(0.0,lalpha);
    }

    //--------------------------------------------------
    //try metrop
    double mul,mur; //means for new bottom nodes, left and right
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      //mul = heterdrawnodemu(bl,Ml,pi.tau,gen);
      //mur = heterdrawnodemu(br,Mr,pi.tau,gen);
      //x.birthp(nx,v,c,mul,mur);
      //x.birthp(nx,v,c,0,0);
      nv[v]++;
      return true;
    } else {
      x = xold;
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal
    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    tree xold = x;
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    x.deathp(nx,0);
    //--------------------------------------------------
    //compute sufficient statistics
    double br,bl; //sums of weights
    double Ml, Mr; //weighted sums of y
    hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma);
    //--------------------------------------------------
    //compute alpha
    //double lhl, lhr, lht;
    //lhl = heterlh(bl,Ml,pi.tau);
    //lhr = heterlh(br,Mr,pi.tau);
    //lht = heterlh(bl+br,Ml+Mr,pi.tau);
    //cout << "flag 1" << r[0] << ' ' << r[di.n-1] <<  endl;
    //x.printtree(xi);
    double lh1 = hetergetmargprob_all(x, xi, di, r, sigma[0], pi, sigma_m, kappa);
    //cout << "flag 2" << endl;
    double lh2 = hetergetmargprob_all(xold, xi, di, r, sigma[0], pi, sigma_m, kappa);
    double lalpha = log(pr) + (lh1 - lh2);
    lalpha = std::min(0.0,lalpha);
    //--------------------------------------------------
    //try metrop
    //double a,b,s2,yb;
    //double mu;
    if(log(gen.uniform()) < lalpha) {
      //mu = heterdrawnodemu(bl+br,Ml+Mr,pi.tau,gen);
      nv[nx->getv()]--;
      //x.deathp(nx,mu);
      //x.deathp(nx,0);
      return true;
    } else {
      x = xold;
      return false;
    }
  }
}



// only for testing: when model is non-spatial
double hetergetmargprob_test(
    tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi
) {
  Rcpp::NumericVector ydata;
  for (int i = 0; i < di.n; i++) {
    ydata.push_back(di.y[i]);
  }
  std::vector<size_t> idv1;
  getbotsid(x, di.x, xi, di.n, di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.n);
  if (df.ncol() > 1) {
    Rcpp::NumericMatrix dfmat = testDFtoNM(df);
    arma::mat Cmat = convertC(dfmat);
    arma::mat temp1 = Cmat.t()*Cmat;
  }
  //for (int i = 0; i < C_.n_cols; i++) {
  //  cout << C_(i,i) << " ";
  //}
  //cout << endl;
  std::set<int> s;
  unsigned size = idv1.size();
  bool flag_true = 1;
  bool flag_false = 0;
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv2;
  idv2.assign( s.begin(), s.end() );
  Rcpp::CharacterVector cnames;
  if (df.size() == 1) {
    df.names() = CharacterVector({"X0"});
    cnames = CharacterVector({"X0"});
  } else {
    cnames = df.names();
  }
  df.push_back(ydata, "y");
  //cout << di.n << " " << di.p << endl;
  Environment env("package:INLA");
  Function inla = env["inla"];
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  std::string formula1 = "y ~ -1 + ";
  for (int i = 0; i < cnames.size(); i++) {
    formula1 += cnames[i];
    if (i != cnames.size()-1) {
      formula1 += " + ";
    }
  }
  //cout << formula1 << endl;
  Function form("as.formula");
  NumericMatrix coords = internal::convert_using_rfunction(dfout, "as.matrix");
  //double logprecstart = -log(pi.sigma2e+sigma*sigma);
  double logprecstart = -2.0*log(sigma);
  // Rcpp::List par_list = Rcpp::List::create(
  //   Rcpp::Named("initial") = logprecstart,
  //   Rcpp::Named("fixed") = flag_true,
  //   Rcpp::Named("prior") = "loggamma",
  //   Rcpp::Named("param") = paramvec
  // );
  Rcpp::List par_list = Rcpp::List::create(
    Rcpp::Named("initial") = logprecstart,
    Rcpp::Named("fixed") = flag_true
  );
  Rcpp::List prec_list = Rcpp::List::create(
    Rcpp::Named("prec") = par_list
  );
  Rcpp::List hyper_list = Rcpp::List::create(
    Rcpp::Named("hyper") = prec_list
  );
  Rcpp::List fixed_list = Rcpp::List::create(
    Rcpp::Named("mean") = 0.0,
    Rcpp::Named("prec") = 1/(pi.tau*pi.tau),
    Rcpp::Named("mean.intercept") = 0.0,
    Rcpp::Named("prec.intercept") = 1/(pi.tau*pi.tau)
  );
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 Rcpp::_["data"] = df,
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["control.fixed"] = fixed_list,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = flag_false,
                                                                  Named("mlik") = flag_true)
                                                                  );
  Rcpp::NumericVector vout = inla_results["mlik"];
  //cout << vout.length();
  //cout << "mlik is: " << vout[0] << endl;
  return vout[0];
  //vout[0];
}

void printX(dinfo& di) {
  double* x = di.x;
  std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt", MyFile.out | MyFile.app);
  for (int i = 0; i < di.n; i++) {
    for (int j = 0; j < di.p; j++) {
      MyFile << *(x + di.p*i+j) << " ";
    }
    MyFile << endl;
  }
  MyFile << "==== spat coords ====" << endl;
  for (int i = 0; i < di.n; i++) {
    MyFile << di.s1[i] << " " << di.s2[i] << endl;
  }
  MyFile.close();
}


// final stage: consider everything
double hetergetmargprob_all(
    tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi,
    double &sigma_m, double &kappa
) {
  //std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt", MyFile.out | MyFile.app);
  Rcpp::NumericVector ydata;
  //MyFile << "==== Data ====" << endl;
  //MyFile << "r[i]:" << endl;
  for (int i = 0; i < di.n; i++) {
    //MyFile << r[i] <<  " " << endl;
    ydata.push_back(r[i]);
  }
  std::vector<size_t> idv1;
  getbotsid(x, di.x, xi, di.n, di.p, idv1);
  //MyFile << "==== X Matrix ====" << endl;
  //printX(di);
  //MyFile << "==== C Matrix ====" << endl;
  Rcpp::DataFrame df = makeC(idv1, di.n);
  //MyFile << endl;
  //printC(df);
  //x.printtreeToFile(xi);
  if (df.ncol() > 1) {
    Rcpp::NumericMatrix dfmat = testDFtoNM(df);
    arma::mat Cmat = convertC(dfmat);
    arma::mat temp1 = Cmat.t()*Cmat;
  }
  std::set<int> s;
  unsigned size = idv1.size();
  bool flag_true = 1;
  bool flag_false = 0;
  // C matrix
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv2;
  idv2.assign( s.begin(), s.end() );
  Rcpp::CharacterVector cnames;
  if (df.size() == 1) {
    df.names() = CharacterVector({"X0"});
    cnames = CharacterVector({"X0"});
    CharacterVector cvec = df.names();
    //Rcpp::NumericVector dv = df[0];
  } else {
    cnames = df.names();
  }
  //df.push_back(ydata, "y");
  //cout << di.n << " " << di.p << endl;
  Environment env("package:INLA");
  Function inla = env["inla"];
  // spatial dataset
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  std::string formula1 = "y ~ -1 + ";
  for (int i = 0; i < cnames.size(); i++) {
    formula1 += cnames[i];
    if (i != cnames.size()-1) {
      formula1 += " + ";
    }
  }
  // adding spatial component
  formula1 += "+f(i, model = spde)";
  //formula1 += "f(i, model = spde)";
  // cout << formula1 << endl;
  //MyFile << "==== Formula ====" << endl;
  //MyFile << formula1 << endl;
  Function form("as.formula");

  NumericMatrix coords = internal::convert_using_rfunction(dfout, "as.matrix");
  //MyFile << "==== spat var, matern var, kappa09 ====" << endl;
  //MyFile << sigma << " " << sigma_m  << " " << kappa << endl;
  //MyFile << sqrt(pi.nu*8.0)/kappa << " " << sigma_m << endl;
  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];

  // SPDE stuff
  Rcpp::NumericVector priorRange = {sqrt(pi.nu*8.0)/kappa,NA_REAL};
  Rcpp::NumericVector priorSigma = {sigma_m,NA_REAL};
  Rcpp::List spde = inlaSpde2PcMatern(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["prior.range"] = priorRange,
    Rcpp::_["prior.sigma"] = priorSigma
  );
  //MyFile << "==== spat var, matern var, kappa10 ====" << endl;
  //MyFile << sigma << " " << sigma_m  << " " << kappa << endl;

  Rcpp::List mesh = di.mesh;
  int nesh = mesh["n"];
  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["loc"] = coords
  );
  //MyFile << "==== spat var, matern var, kappa11 ====" << endl;
  //MyFile << sigma << " " << sigma_m  << " " << kappa << endl;

  Function inlaStack = env["inla.stack"];
  Rcpp::List y_list = Rcpp::List::create(Named("y") = ydata);
  Rcpp::List A_list = Rcpp::List::create(
    Named("A") = A,
    Named("intercept") = 1);
  int nspde = spde["n.spde"];
  //MyFile << "nesh: " << nesh <<  "nspde: " <<  nspde << endl;
  Rcpp::List effect_list;
  if (df.size() == 1) {
    //MyFile << "onevar" << endl;
    effect_list = Rcpp::List::create(
      //Rcpp::List::create(Named("i") = Rcpp::seq(1,nspde)),
      Named("i") = Rcpp::seq(1,nspde),
      Named("X0") = df[0]//df[0]
    );
  } else {
    effect_list = Rcpp::List::create(
      //Rcpp::List::create(Named("i") = Rcpp::seq(1,nspde)),
      Named("i") = Rcpp::seq(1,nspde),
      df
    );
  }
  Rcpp::List stkDat = inlaStack(
    Rcpp::_["data"] = y_list,
    Rcpp::_["A"] = A_list,
    Rcpp::_["effects"] = effect_list
  );
  Function inlaStackData = env["inla.stack.data"];
  Function inlaStackA = env["inla.stack.A"];
  Rcpp::List controlpred_list = Rcpp::List::create(
    Rcpp::Named("A") = inlaStackA(stkDat),
    Rcpp::Named("compute") = flag_true
  );
  // residual variance part
  double logprecstart = -2.0*log(sigma);
  Rcpp::NumericVector paramvec ={2.5,0.25};
  Rcpp::List par_list = Rcpp::List::create(
     Rcpp::Named("initial") = logprecstart,
     Rcpp::Named("fixed") = flag_true,
     Rcpp::Named("prior") = "loggamma",
     Rcpp::Named("param") = paramvec
   );
  Rcpp::List prec_list = Rcpp::List::create(
    Rcpp::Named("prec") = par_list
  );
  Rcpp::List hyper_list = Rcpp::List::create(
    Rcpp::Named("hyper") = prec_list
  );
  Rcpp::List fixed_list = Rcpp::List::create(
    Rcpp::Named("mean") = 0.0,
    Rcpp::Named("prec") = 1/(pi.tau*pi.tau),
    Rcpp::Named("mean.intercept") = 0.0,
    Rcpp::Named("prec.intercept") = 1/(pi.tau*pi.tau)
  );
  //cout << "pitau" << pi.tau << endl;
  //MyFile << "==== spat var, matern var, kappa ====" << endl;
  //MyFile << "sigma: " << sigma << endl;
  //MyFile << sigma << " " << sigma_m  << " " << kappa << endl;
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["control.fixed"] = fixed_list,
                                 Rcpp::_["num.threads"] = "1:1",
                                 Rcpp::_["only.hyperparam"] = flag_true,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = flag_false,
                                                                  Named("mlik") = flag_true),
                                 Rcpp::_["data"] = inlaStackData(stkDat, Rcpp::_["spde"] = spde),
                                 Rcpp::_["control.predictor"] = controlpred_list);
  //MyFile << "==== marglik ====" << endl;
  Rcpp::NumericVector vout = inla_results["mlik"];
  //cout << vout.length();
  //cout << "mlik is: " << vout[0] << endl;
  //MyFile << vout[0] << endl;
   return vout[0];
  //MyFile.close();
  //vout[0];
}



// for test only
double hetergetmargprob_testonly(
    tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi
) {
  Rcpp::NumericVector ydata;
  for (int i = 0; i < di.n; i++) {
    ydata.push_back(r[i]);
  }
  std::vector<size_t> idv1;
  getbotsid(x, di.x, xi, di.n, di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.n);
  Rcpp::NumericMatrix dfmat = testDFtoNM(df);
  arma::mat Cmat = convertC(dfmat);
  arma::mat temp1 = Cmat*Cmat.t()*pi.tau*pi.tau;
  arma::mat temp2(Cmat.n_rows , Cmat.n_rows, fill::eye);
  arma::vec ydat(ydata);
  arma::mat temp3 = temp2 * sigma * sigma;
  arma::mat covs = temp1 + temp3;
  //cout << "determinant is " << arma::log_det_sympd(covs) << endl;
  double lh = dmvnorm_rcpp(ydat , covs);
  //cout << lh << endl;
  return lh;
}


double hetergetmargprob(
    tree& x, xinfo& xi, dinfo& di, double* r, double sigma, pinfo& pi
) {
  Rcpp::NumericVector ydata;
  for (int i = 0; i < di.n; i++) {
    ydata.push_back(r[i]);
  }
  std::vector<size_t> idv1;
  getbotsid(x, di.x, xi, di.n, di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.n);
  std::set<int> s;
  unsigned size = idv1.size();
  bool flag_true = 1;
  bool flag_false = 0;
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv2;
  idv2.assign( s.begin(), s.end() );
  Rcpp::CharacterVector cnames;
  if (df.size() == 1) {
    df.names() = CharacterVector({"X0"});
    cnames = CharacterVector({"X0"});
  } else {
    cnames = df.names();
  }
  //cout << di.n << " " << di.p << endl;
  Environment env("package:INLA");
  Function inla = env["inla"];

  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);

  std::string formula1 = "y ~ -1 + ";
  for (int i = 0; i < cnames.size(); i++) {
    formula1 += cnames[i];
    formula1 += " + ";
  }
  formula1 += "f(i, model = spde)";
  cout << formula1 << endl;
  Function form("as.formula");
  /*
  Function inlaMesh2D = env["inla.mesh.2d"];
  Rcpp::NumericVector edge_arg = {0.1,0.3};
  Rcpp::NumericVector cutoff = {0.06};

  Rcpp::List mesh = inlaMesh2D(Rcpp::_["loc"] = dfout,
                               Rcpp::_["max.edge"] = edge_arg,
                               Rcpp::_["cutoff"] = cutoff);
  */

  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];
  Rcpp::NumericVector priorRange = {sqrt(pi.nu*8)/pi.kappa,NA_REAL};
  Rcpp::NumericVector priorSigma = {sqrt(pi.sigma2x),NA_REAL};
  Rcpp::List spde = inlaSpde2PcMatern(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["prior.range"] = priorRange,
    Rcpp::_["prior.sigma"] = priorSigma
  );

  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  NumericMatrix coords = internal::convert_using_rfunction(dfout, "as.matrix");
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["loc"] = coords
  );

  Function inlaStack = env["inla.stack"];
  Rcpp::List y_list = Rcpp::List::create(Named("y") = ydata);
  Rcpp::List A_list = Rcpp::List::create(
    Named("A") = A,
    Named("intercept") = 1);
  int nspde = spde["n.spde"];
  Rcpp::List effect_list = Rcpp::List::create(
    Rcpp::List::create(Named("i") = Rcpp::seq(1,nspde)),
    Rcpp::List::create(Named("Intercept") = 1,
                       Named("coords") = df)
  );
  Rcpp::List stkDat = inlaStack(
    Rcpp::_["data"] = y_list,
    Rcpp::_["A"] = A_list,
    Rcpp::_["effects"] = effect_list
  );
  Function inlaStackData = env["inla.stack.data"];
  Function inlaStackA = env["inla.stack.A"];
  Rcpp::List controlpred_list = Rcpp::List::create(
    Rcpp::Named("A") = inlaStackA(stkDat),
    Rcpp::Named("compute") = flag_true
  );

  double logprecstart = -2.0*log(pi.sigma2e+0.05*0.05);

  Rcpp::NumericVector paramvec ={1.0,0.1};
  Rcpp::List par_list = Rcpp::List::create(
    Rcpp::Named("initial") = logprecstart,
    Rcpp::Named("fixed") = flag_true,
    Rcpp::Named("prior") = "loggamma",
    Rcpp::Named("param") = paramvec
  );

  Rcpp::List prec_list = Rcpp::List::create(
    Rcpp::Named("prec") = par_list
  );

  Rcpp::List hyper_list = Rcpp::List::create(
    Rcpp::Named("hyper") = prec_list
  );
 Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = flag_false,
                                                                  Named("mlik") = flag_true),
                                 Rcpp::_["data"] = inlaStackData(stkDat, Rcpp::_["spde"] = spde),
                               Rcpp::_["control.predictor"] = controlpred_list);
  Rcpp::NumericVector vout = inla_results["mlik"];
  //cout << vout.length();
  cout << "mlik is: " << vout[0] << endl;
  return vout[0];
  //vout[0];
}



bool heterbd_new(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double sigma,
             std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
  if(gen.uniform() < PBx) { //do birth or death
    //--------------------------------------------------
    //draw proposal
    cout << "birth" << endl;
    tree::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    tree xold = x;
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    size_t nr,nl; //counts in proposed bots
    countbotnodes(x, nx, v, c, xi, di, nl, nr);
    cout << "new tree"<< endl;
    x.birthp(nx,v,c,0,0);
    x.printtree(xi);
    cout << "old tree"<< endl;
    xold.printtree(xi);
    //--------------------------------------------------
    //compute sufficient statistics
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0;
    cout << nl << " " << nr << endl;
    if((nl>=5) && (nr>=5)) { //cludge?
      double lh1, lh2;
      lh1 = hetergetmargprob(x, xi, di, r, sigma, pi);
      cout << "flag 1 " << endl;
      lh2 = hetergetmargprob(xold, xi, di, r, sigma, pi);
      cout << "flag 2 " << endl;

      alpha=1.0;
      lalpha = log(pr) + (lh1-lh2);
      lalpha = std::min(0.0,lalpha);
    }
    //--------------------------------------------------
    //try metrop
    //double mul,mur; //means for new bottom nodes, left and right
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      cout << "new tree" << endl;
      nv[v]++;
      cout << "final tree"<< endl;
      x.printtree(xi);
      return true;
    } else {
      x = xold;
      cout << "final tree"<< endl;
      x.printtree(xi);
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal
    cout << "death" << endl;
    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    tree xold = x;
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    x.deathp(nx,0);
    x.printtree(xi);
    xold.printtree(xi);
    double lh1, lh2;
    lh1 = hetergetmargprob(x, xi, di, r, sigma, pi);
    lh2 = hetergetmargprob(xold, xi, di, r, sigma, pi);
    double lalpha = log(pr) + (lh1-lh2);
    lalpha = std::min(0.0,lalpha);
    //--------------------------------------------------
    //try metrop
    if(log(gen.uniform()) < lalpha) {
      cout << "new tree" << endl;
      nv[nx->getv()]--;
      return true;
    } else {
      x = xold;
      return false;
    }
  }
}


bool heterbd_new_test(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double sigma,
                 std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
  if(gen.uniform() < PBx) { //do birth or death
    //--------------------------------------------------
    //draw proposal
    //cout << "birth" << endl;
    tree::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    tree xold = x;
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    size_t nr,nl; //counts in proposed bots
    double Ml, Mr, bl, br;
    double sigmas[di.n];
    for (int i=0; i <di.n; i++) {
      sigmas[i] = sigma;
    }
    hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigmas,r);
    //countbotnodes(x, nx, v, c, xi, di, nl, nr);
    //cout << "new tree"<< endl;
    //cout << "sigma: " << sigma << " tau: " << pi.tau << endl;
    x.birthp(nx,v,c,0,0);
    //x.printtree(xi);
    //cout << "old tree"<< endl;
    //xold.printtree(xi);
    //--------------------------------------------------
    //compute sufficient statistics
    //hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma);
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0;
    //cout << nl << " " << nr << endl;
    if((nl>=5) && (nr>=5)) { //cludge?
      double lh1, lh2, lhl, lhr, lht;
      lh1 = hetergetmargprob_test(x, xi, di, r, sigma, pi);
      cout << "flag1" << endl;
      lh2 = hetergetmargprob_test(xold, xi, di, r, sigma, pi);
      cout << "flag2" << endl;
      //lhl = heterlh(bl,Ml,pi.tau);
      //lhr = heterlh(br,Mr,pi.tau);
      //lht = heterlh(bl+br,Ml+Mr,pi.tau);

      alpha=1.0;
      lalpha = log(pr) + lh1 - lh2;
      //lalpha = log(pr) + (lhl + lhr - lht);
      //cout << "delta: " << (lh1 - lh2) - (lhl + lhr - lht) << endl;
      //cout << "loglh: " << (lh1 - lh2);
      //cout << "loglh2: " << (lhl + lhr - lht);
      lalpha = std::min(0.0,lalpha);
    }
    //--------------------------------------------------
    //try metrop
    //double mul,mur; //means for new bottom nodes, left and right
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      //cout << "new tree" << endl;
      nv[v]++;
      //cout << "final tree"<< endl;
      //x.printtree(xi);
      return true;
    } else {
      x = xold;
      //cout << "final tree"<< endl;
      //x.printtree(xi);
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal
    //cout << "death" << endl;
    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    tree xold = x;
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    double bl, br, Ml, Mr;
    size_t nl, nr;
    double sigmas[di.n];
    for (int i=0; i <di.n; i++) {
      sigmas[i] = sigma;
    }
    hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigmas);
    x.deathp(nx,0);
    //x.printtree(xi);
    //xold.printtree(xi);
    double lh1, lh2, lhl, lhr, lht;
    lh1 = hetergetmargprob_test(x, xi, di, r, sigma, pi);
    lh2 = hetergetmargprob_test(xold, xi, di, r, sigma, pi);
    //lhl = heterlh(bl,Ml,pi.tau);
    //lhr = heterlh(br,Mr,pi.tau);
    //lht = heterlh(bl+br,Ml+Mr,pi.tau);

    //alpha=1.0;
    //cout << "delta: " << (lh1 - lh2) + (lhl + lhr - lht) << endl;
    //cout << "loglh: " << (lh1 - lh2);
    //cout << "loglh2: " << -(lhl + lhr - lht);
    double lalpha = log(pr) + lh1 - lh2;
    //double lalpha = log(pr) -(lhl + lhr - lht);
    lalpha = std::min(0.0,lalpha);
    //--------------------------------------------------
    //try metrop
    if(log(gen.uniform()) < lalpha) {
      //cout << "new death tree" << endl;
      nv[nx->getv()]--;
      //x.printtree(xi);
      return true;
    } else {
      x = xold;
      //x.printtree(xi);
      return false;
    }
  }
}

bool heterbd_new_test2(tree& x, xinfo& xi, double* r, dinfo& di, pinfo& pi, double sigma,
                      std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
  if(gen.uniform() < PBx) { //do birth or death
    //--------------------------------------------------
    //draw proposal
    //cout << "birth" << endl;
    tree::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    tree xold = x; // a new tree
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    size_t nr,nl; //counts in proposed bots
    double Ml, Mr, bl, br;
    double sigmas[di.n];
    for (int i=0; i <di.n; i++) {
      sigmas[i] = sigma;
    }
    hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigmas,r);
    //countbotnodes(x, nx, v, c, xi, di, nl, nr);
    //cout << "new tree"<< endl;
    //cout << "sigma: " << sigma << " tau: " << pi.tau << endl;
    x.birthp(nx,v,c,0,0);
    //x.printtree(xi);
    //cout << "old tree"<< endl;
    //xold.printtree(xi);
    //--------------------------------------------------
    //compute sufficient statistics
    //hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma);
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0;
    //cout << nl << " " << nr << endl;
    if((nl>=5) && (nr>=5)) { //cludge?
      double lh1, lh2, lhl, lhr, lht;
      //lh1 = hetergetmargprob_test(x, xi, di, r, sigma, pi);
      //lh2 = hetergetmargprob_test(xold, xi, di, r, sigma, pi);
      lhl = heterlh(bl,Ml,pi.tau);
      lhr = heterlh(br,Mr,pi.tau);
      lht = heterlh(bl+br,Ml+Mr,pi.tau);
      alpha=1.0;
      lalpha = log(pr) + (lhl + lhr - lht);
      //cout << "delta: " << (lh1 - lh2) - (lhl + lhr - lht) << endl;
      //cout << "loglh: " << (lh1 - lh2);
      //cout << "loglh2: " << (lhl + lhr - lht);
      lalpha = std::min(0.0,lalpha);
    }
    //--------------------------------------------------
    //try metrop
    //double mul,mur; //means for new bottom nodes, left and right
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      //x.birthp(nx,v,c,0,0);
      //x = xnew;
      //x = xnew;
      //cout << "new tree" << endl;
      nv[v]++;
      //cout << "final tree"<< endl;
      //x.printtree(xi);
      return true;
    } else {
      x = xold;
      //cout << "final tree"<< endl;
      //x.printtree(xi);
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal
    //cout << "death" << endl;
    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    tree xold = x;
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    double bl, br, Ml, Mr;
    size_t nl, nr;
    double sigmas[di.n];
    for (int i=0; i <di.n; i++) {
      sigmas[i] = sigma;
    }
    hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigmas);
    x.deathp(nx,0);
    //x.printtree(xi);
    //xold.printtree(xi);
    double lh1, lh2, lhl, lhr, lht;
    //lh1 = hetergetmargprob_test(x, xi, di, r, sigma, pi);
    //lh2 = hetergetmargprob_test(xold, xi, di, r, sigma, pi);
    lhl = heterlh(bl,Ml,pi.tau);
    lhr = heterlh(br,Mr,pi.tau);
    lht = heterlh(bl+br,Ml+Mr,pi.tau);

    //alpha=1.0;
    //cout << "delta: " << (lh1 - lh2) + (lhl + lhr - lht) << endl;
    //cout << "loglh: " << (lh1 - lh2);
    //cout << "loglh2: " << -(lhl + lhr - lht);
    double lalpha = log(pr) -(lhl + lhr - lht);
    lalpha = std::min(0.0,lalpha);
    //--------------------------------------------------
    //try metrop
    if(log(gen.uniform()) < lalpha) {
      //cout << "new death tree" << endl;
      //x.deathp(nx,0);
      //x = xold;
      nv[nx->getv()]--;
      //x.printtree(xi);
      return true;
    } else {
      x = xold;
      //x.printtree(xi);
      return false;
    }
  }
}

double heterbd_drawsigma_margprob(double* r, dinfo& di, pinfo& pi, double sigma, double &sigma_m, double &kappa) {
  Rcpp::NumericVector ydata;
  //cout << "---data---" << endl;
  for (int i = 0; i < di.n; i++) {
    ydata.push_back(r[i]);
  //  cout << r[i] << endl;
  }
  Environment env("package:INLA");
  Function inla = env["inla"];
  bool flag_true = 1;
  bool flag_false = 0;
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  std::string formula1 = "y ~ -1 + ";
  formula1 += "f(i, model = spde)";
  //cout << formula1 << endl;
  Function form("as.formula");
  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];
  //cout << kappa << "aaa " << sigma_m << endl;
  Rcpp::NumericVector priorRange = {sqrt(pi.nu*8.0)/kappa,NA_REAL};
  Rcpp::NumericVector priorSigma = {sigma_m,NA_REAL};
  Rcpp::List spde = inlaSpde2PcMatern(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["prior.range"] = priorRange,
    Rcpp::_["prior.sigma"] = priorSigma
  );
  Rcpp::List mesh = di.mesh;
  int nesh = mesh["n"];
  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  NumericMatrix coords = internal::convert_using_rfunction(dfout, "as.matrix");
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["loc"] = coords
  );
  Function inlaStack = env["inla.stack"];
  Rcpp::List y_list = Rcpp::List::create(Named("y") = ydata);
  Rcpp::List A_list = Rcpp::List::create(
    Named("A") = A,
    Named("intercept") = 1);
  int nspde = spde["n.spde"];
  //cout << nspde << endl;
  Rcpp::List effect_list = Rcpp::List::create(
    Rcpp::List::create(Named("i") = Rcpp::seq(1,nspde)),
    Rcpp::List::create(Named("Intercept") = Rcpp::rep(1, di.n))
  );
  Rcpp::List stkDat = inlaStack(
    Rcpp::_["data"] = y_list,
    Rcpp::_["A"] = A_list,
    Rcpp::_["effects"] = effect_list
  );
  Function inlaStackData = env["inla.stack.data"];
  Function inlaStackA = env["inla.stack.A"];
  Rcpp::List controlpred_list = Rcpp::List::create(
    Rcpp::Named("A") = inlaStackA(stkDat),
    Rcpp::Named("compute") = flag_true
  );
  double logprecstart = -2.0*log(sigma);
  Rcpp::NumericVector paramvec ={2.5,0.25};
  Rcpp::List par_list = Rcpp::List::create(
    Rcpp::Named("initial") = logprecstart,
    Rcpp::Named("fixed") = flag_true,
    Rcpp::Named("prior") = "loggamma",
    Rcpp::Named("param") = paramvec
  );
  Rcpp::List prec_list = Rcpp::List::create(
    Rcpp::Named("prec") = par_list
  );

  Rcpp::List hyper_list = Rcpp::List::create(
    Rcpp::Named("hyper") = prec_list
  );
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 //Rcpp::_["scale"] = di.dat_size,
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["num.threads"] = "1:1",
                                 Rcpp::_["only.hyperparam"] = flag_true,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = flag_false,
                                                                  Named("mlik") = flag_true),
                                                                  Rcpp::_["data"] = inlaStackData(stkDat, Rcpp::_["spde"] = spde),
                                                                  Rcpp::_["control.predictor"] = controlpred_list);
  Rcpp::NumericVector vout = inla_results["mlik"];
  //cout << sigma << " " << sigma_m << " " << kappa << " " << di.n << endl;
  //cout << vout[0] << endl;
  //cout << vout.length();
  //cout << "mlik is: " << vout[0] << endl;
  return vout[0];
}

// testing for Gaussian code

double heterbd_drawsigma_margprob_test(double* r, dinfo& di, pinfo& pi, double sigma) {
  Rcpp::NumericVector ydata;
  for (int i = 0; i < di.n; i++) {
    ydata.push_back(r[i]);
  }
  DataFrame df = DataFrame::create(Named("y") = ydata);
  Environment env("package:INLA");
  Function inla = env["inla"];
  bool flag_true = 1;
  bool flag_false = 0;
  Rcpp::CharacterVector namevec;
  std::string formula1 = "y ~ -1";
  //cout << formula1 << endl;
  Function form("as.formula");
  //Rcpp::NumericVector priorRange = {sqrt(pi.nu*8)/pi.kappa,NA_REAL};
  //Rcpp::NumericVector priorSigma = {sqrt(pi.sigma2x),NA_REAL};
  Rcpp::List y_list = Rcpp::List::create(Named("y") = ydata);
  double logprecstart = -2.0*log(sigma);
  Rcpp::List par_list = Rcpp::List::create(
    Rcpp::Named("initial") = logprecstart,
    Rcpp::Named("fixed") = flag_true
  );
  Rcpp::List prec_list = Rcpp::List::create(
    Rcpp::Named("prec") = par_list
  );
  Rcpp::List hyper_list = Rcpp::List::create(
    Rcpp::Named("hyper") = prec_list
  );
  Rcpp::List fixed_list = Rcpp::List::create(
    Rcpp::Named("mean") = 0.0,
    Rcpp::Named("prec") = 1/(pi.tau*pi.tau),
    Rcpp::Named("mean.intercept") = 0.0,
    Rcpp::Named("prec.intercept") = 1/(pi.tau*pi.tau)
  );
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 Rcpp::_["data"] = df,
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["control.fixed"] = fixed_list,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = flag_false,
                                                                  Named("mlik") = flag_true)
  );
  Rcpp::NumericVector vout = inla_results["mlik"];
  //cout << vout.length();
  //cout << "mlik is: " << vout[0] << endl;
  //cout << "first residual: " << r[0] << endl;
  //double trueloglik = 0.0;
  //for (int i = 0; i < di.n; i++) {
  //  trueloglik -= r[i]*r[i]/2.0/sigma/sigma;
  //}
  //trueloglik -= di.n/2.0*log(sigma*sigma*2*PI);
  //cout << "numer mlik is: " << trueloglik << endl;
  return vout[0];
}


double heterbd_drawsigma(
    dinfo& di, double* r, double sigma,double &kappa,double &sigma_m, pinfo& pi, double lambda, double nu
) {
  arma::arma_rng::set_seed(1);
  // proposing
  // check if this is doing what we intended
  double omi = 0.08; // propose parameter for MH draws
  double rands = arma::randn<double>();
  double temp = rands*omi+(2.0)*std::log(sigma);
  double sigma_star = std::sqrt(std::exp(temp));
  //double sigma_star = rands*omi + sigma;
  //if (sigma_star < 0) {
  //  return sigma;
  //} else {
  //cout << "sigma, sigmam, kappa: " << sigma_star << " " << sigma_m << " " << kappa << endl;
  double mlik1 = heterbd_drawsigma_margprob(r, di, pi, sigma_star, sigma_m, kappa);
  //cout << "sigma, sigmam, kappa: " << sigma << " " << sigma_m << " " << kappa << endl;
  double mlik2 = heterbd_drawsigma_margprob(r, di, pi, sigma, sigma_m, kappa);
  bool flag_true = 1;
  //cout << lambda << "lambda" << nu << "nu" << endl;
  double p1 = R::dchisq(nu*lambda/sigma_star/sigma_star,nu,flag_true);
  //cout << p1 << "p1" << endl;
  double p2 = R::dchisq(nu*lambda/sigma/sigma,nu,flag_true);
  double pdiff = -nu*lambda/(2.0)*(1/sigma_star/sigma_star - 1/sigma/sigma) + (1.0 + nu/2.0)*(std::log(sigma) - std::log(sigma_star))*(2.0);
  double pdiff2 = (p1 - p2) + (std::log(sigma) - std::log(sigma_star))*(4.0);
  //cout << pdiff << " " << pdiff2 << endl;
  //cout << p1 << " " <<  nu*lambda/sigma_star/sigma_star << " " << nu << endl;
  //cout << sigma_star-sigma << " sigma " <<  mlik1- mlik2 << " mlik " << endl;
  double mh_ratio = mlik1 - mlik2 + (p1 - p2) + (std::log(sigma_star) - std::log(sigma))*(2.0);
  //double mh_ratio = std::min(0.0, mlik1 - mlik2  -(p1 - p2));
  //cout << mh_ratio << " " << sigma_star << " " << sigma << endl;
  if (mh_ratio > 0) {
    return sigma_star;
   } else {
    double random_gen = arma::randu<double>();
    //cout << random_gen << endl;
    if (random_gen < std::exp(mh_ratio)) {
      return sigma_star;
    }
    return sigma;
  }
  //}
}

void heterbd_drawspathyperpars(
    dinfo& di, double* r, double sigma, double &kappa, double &sigma_m, pinfo& pi, double &mlik) {
  arma::arma_rng::set_seed(1);
  // proposing
  // check if this is doing what we intended
  double omi = 0.08; // propose parameter for MH draws
  double rands = arma::randn<double>();
  double rands2 = arma::randn<double>();
  double temp_kappa = rands*omi+std::log(kappa);
  double temp_sigma_m = rands2*omi+std::log(sigma_m);
  double sigma_mstar = std::exp(temp_sigma_m);
  double kappa_star = std::exp(temp_kappa);
  //cout << "sigma, sigmam, kappa: " << sigma << " " << sigma_mstar << " " << kappa_star << endl;
  double mlik1 = heterbd_drawsigma_margprob(r, di, pi, sigma, sigma_mstar, kappa_star);
  //cout << "sigma, sigmam, kappa: " << sigma << " " << sigma_m << " " << kappa << endl;
  double mlik2 = heterbd_drawsigma_margprob(r, di, pi, sigma, sigma_m, kappa);
  //cout << mlik2;
  // pc prior difference
  //cout << pi.sigma_m0 << " " << pi.alpha_2 << endl;
  double rho = 2.0*std::sqrt(2.0*pi.nu)/kappa;
  double rho_star = 2.0*std::sqrt(2.0*pi.nu)/kappa_star;
  double pdiff = std::log(pi.alpha_1)*pi.rho_0*(1/rho_star - 1/rho) + std::log(pi.alpha_2)/pi.sigma_m0*(sigma_mstar - sigma_m) - (2.0)*(std::log(rho_star) - std::log(rho));
  bool flag_true = 1;
  double mh_ratio = mlik1 - mlik2 + pdiff + (std::log(sigma_mstar) - std::log(sigma_m)) + (std::log(kappa_star) - std::log(kappa));
  if (mh_ratio > 0) {
    kappa = kappa_star;
    sigma_m = sigma_mstar;
  } else {
    double random_gen = arma::randu<double>();
    if (random_gen < std::exp(mh_ratio)) {
      kappa = kappa_star;
      sigma_m = sigma_mstar;
    }
  }
}


// weighted INLA using agaussian
double mlikcalcweighted (double* r, dinfo& di, pinfo& pi, double sigma, double &sigma_m, double &kappa) {
  Rcpp::NumericVector r_;
  for (int i = 0; i < di.n; i++) {
    r_.push_back(r[i]);
  }
  Rcpp::NumericVector size_;
  for (int i = 0; i < di.n; i++) {
    size_.push_back(di.dat_size[i]);
  }
  arma::vec rvec(r_);
  arma::vec sizevec(size_);
  Environment env("package:INLA");
  Environment bas("package:base");
  Function inla = env["inla"];
  Function inlaag = env["inla.agaussian"];
  Function unclass = bas["unclass"];
  r_.attr("dim") = Dimension(di.n, 1);
  size_.attr("dim") = Dimension(di.n, 1);
  Rcpp::NumericMatrix mr_ = as<NumericMatrix>(r_);
  Rcpp::NumericMatrix msize_ = as<NumericMatrix>(size_);
  //cout << mr_.rows() << endl;
  Rcpp::List Y = inlaag(
    Rcpp::_["y"] = r_,
    Rcpp::_["s"] = msize_);
  Rcpp::List unclassY = unclass(Y);
  Rcpp::NumericVector _Y1;
  Rcpp::NumericVector _Y2;
  Rcpp::NumericVector _Y3;
  Rcpp::NumericVector _Y4;
  Rcpp::NumericVector _Y5;
  Rcpp::NumericVector V1 = unclassY["Y1"];
  Rcpp::NumericVector V2 = unclassY["Y2"];
  Rcpp::NumericVector V3 = unclassY["Y3"];
  Rcpp::NumericVector V4 = unclassY["Y4"];
  Rcpp::NumericVector V5 = unclassY["Y5"];
  for (int i = 0; i < di.n; i++) {
    _Y1.push_back(0.0);
    _Y2.push_back(V2[i]);
    _Y3.push_back(V3[i]);
    _Y4.push_back(V4[i]);
    _Y5.push_back(V5[i]);
  }
  Rcpp::DataFrame Ydf = Rcpp::DataFrame::create(Named("Y1") = _Y1,
                                                Named("Y2") = _Y2,
                                                Named("Y3") = _Y3,
                                                Named("Y4") = _Y4,
                                                Named("Y5") = _Y5);
  bool flag_true = 1;
  bool flag_false = 0;
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  std::string formula1 = "Y ~ -1 + ";
  formula1 += "f(i, model = spde)";
  //cout << formula1 << endl;
  Function form("as.formula");
  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];
  Rcpp::NumericVector priorRange = {sqrt(pi.nu*8.0)/kappa,NA_REAL};
  Rcpp::NumericVector priorSigma = {sigma_m,NA_REAL};
  Rcpp::List spde = inlaSpde2PcMatern(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["prior.range"] = priorRange,
    Rcpp::_["prior.sigma"] = priorSigma
  );
  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  NumericMatrix coords = internal::convert_using_rfunction(dfout, "as.matrix");
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["loc"] = coords
  );
  Function inlaStack = env["inla.stack"];
  Rcpp::List y_list = Rcpp::List::create(Named("Y") = Ydf);
  Rcpp::List A_list = Rcpp::List::create(
    Named("A") = A//,
    //Named("intercept") = 1
  );
  int nspde = spde["n.spde"];
  Rcpp::List effect_list = Rcpp::List::create(
    Named("i") = Rcpp::seq(1,nspde)//,
    //Rcpp::List::create(Named("Intercept") = Rcpp::rep(1, di.n))
  );
  Rcpp::List stkDat = inlaStack(
    Rcpp::_["data"] = y_list,
    Rcpp::_["A"] = A_list,
    Rcpp::_["effects"] = effect_list
  );
  Function inlaStackData = env["inla.stack.data"];
  Function inlaStackA = env["inla.stack.A"];
  Rcpp::List controlpred_list = Rcpp::List::create(
    Rcpp::Named("A") = inlaStackA(stkDat),
    Rcpp::Named("compute") = flag_true
  );
  double logprecstart = -2.0*log(sigma);
  Rcpp::NumericVector paramvec ={1.0,0.1};
  Rcpp::List par_list = Rcpp::List::create(
    Rcpp::Named("initial") = logprecstart,
    Rcpp::Named("fixed") = flag_true,
    Rcpp::Named("prior") = "loggamma",
    Rcpp::Named("param") = paramvec
  );
  Rcpp::List prec_list = Rcpp::List::create(
    Rcpp::Named("prec") = par_list
  );
  Rcpp::List hyper_list = Rcpp::List::create(
    Rcpp::Named("hyper") = prec_list
  );
  Rcpp::List iData = inlaStackData(stkDat, Rcpp::_["spde"] = spde);
  iData["Y"] = Y;
  Rcpp::List tt = iData["Y"];
  Rcpp::NumericVector tp = tt["Y1"];
  // due to numerical issues we need to manually replace the vector as 0
  for (int i = 0; i < tp.length(); i++) {
    tp[i] = 0.0;
    //cout << tp[i] << " ";
  }
  tt["Y1"] = tp;
  iData["Y"] = tt;
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 Rcpp::_["family"] = "agaussian",
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = flag_false,
                                                                  Named("mlik") = flag_true),
                                                                  Rcpp::_["data"] = iData,
                                                                  Rcpp::_["control.predictor"] = controlpred_list);
  Rcpp::NumericVector vout = inla_results["mlik"];
  cout << sigma << " " << sigma_m << " " << kappa << endl;
  cout << vout[0] << endl;
  //Rcpp::DataFrame Ydf = asdf(Yalg);
  //return 0;
  return vout[0];
}

// trivariate hyperparameter update for model checking
void heterbd_trihyperpars(
    dinfo& di, double* r, double& sigma, double &kappa, double &sigma_m, pinfo& pi, double &mlik, double lambda, double nu
) {
  // storing arma vector
  std::vector<double> r_;
  for (int i = 0; i < di.n; i++) {
    r_.push_back(r[i]);
   ///cout << r[i] << endl;
  }
  arma::vec rvec(r_);
  // random walk proposals
  double omi = 0.1; // propose parameter for MH draws
  double rands = arma::randn<double>();
  double rands2 = arma::randn<double>();
  double rands3 = arma::randn<double>();
  double temp_kappa = rands*omi+std::log(kappa);
  double temp_sigma_m = rands2*omi+(2.0)*std::log(sigma_m);
  double temp = rands3*omi+(2.0)*std::log(sigma);
  // propose new spatial std
  double sigma_mstar = std::sqrt(std::exp(temp_sigma_m));
  // propose new Matern scale
  double kappa_star = std::exp(temp_kappa);
  // propose new nonspatial std
  double sigma_star = std::sqrt(std::exp(temp));
  // prior difference
  // spatial part
  double rho = 2.0*std::sqrt(2.0*pi.nu)/kappa;
  double rho_star = 2.0*std::sqrt(2.0*pi.nu)/kappa_star;
  double pdiff = std::log(pi.alpha_1)*pi.rho_0*(1/rho_star - 1/rho) + std::log(pi.alpha_2)/pi.sigma_m0*(sigma_mstar - sigma_m) + (2.0)*(std::log(kappa) - std::log(kappa_star));
  // nonspatial part
  bool flag_true = 1;
  double p1 = R::dchisq(nu*lambda/sigma_star/sigma_star,nu,flag_true);
  double p2 = R::dchisq(nu*lambda/sigma/sigma,nu,flag_true);
  double pdiff_nonspat = p1 - p2;
  // Jacobian parts
  double jacob_spat = (std::log(sigma_mstar) - std::log(sigma_m))*(2.0) + (std::log(kappa_star) - std::log(kappa));
  double jacob_nonspat = (std::log(sigma_star) - std::log(sigma))*(2.0);
  // likelihood parts
  arma::mat sig = calcMat_weighted(di, pi, sigma, sigma_m, kappa);
  arma::mat sig_star = calcMat_weighted(di, pi, sigma_star, sigma_mstar, kappa_star);
  // calculate likelihood
  //Environment env("base");
  //Function determinant = env["determinant"];
  //Rcpp::List temps1 = determinant(Rcpp::_["x"] = sig,
  //                            Rcpp::_["logarithm"] = flag_true);
  //Rcpp::List temps2 = determinant(Rcpp::_["x"] = sig_star,
  //                            Rcpp::_["logarithm"] = flag_true);
  //Rcpp::NumericVector tempss1 = temps1["modulus"];
  //Rcpp::NumericVector tempss2 = temps2["modulus"];
  //double tempsss1 = tempss1[0];
  //double tempsss2 = tempss2[0];
  double mlik1 = -0.5*arma::log_det_sympd(sig) - 0.5*arma::sum(rvec.t()*sig.i()*rvec);
  double mlik2 = -0.5*arma::log_det_sympd(sig_star) - 0.5*arma::sum(rvec.t()*sig_star.i()*rvec);
  double mh_ratio = mlik2 - mlik1 + jacob_spat + jacob_nonspat + pdiff + pdiff_nonspat;
  if (mh_ratio > 0) {
    kappa = kappa_star;
    sigma_m = sigma_mstar;
    sigma = sigma_star;
  } else {
    double random_gen = arma::randu<double>();
    if (random_gen < std::exp(mh_ratio)) {
      kappa = kappa_star;
      sigma_m = sigma_mstar;
      sigma = sigma_star;
    }
  }
}



double dmvnorm_distance( arma::vec y, arma::mat Sigma )
{
  int n = Sigma.n_rows;
  double res=0;
  double fac=1;
  for (int ii=0; ii<n; ii++){
    for (int jj=ii; jj<n; jj++){
      if (ii==jj){ fac = 1; } else { fac = 2;}
      res += fac *y[ii] * Sigma(ii,jj) * y[jj];
    }
  }
  return res;
}

// I used code
double dmvnorm_rcpp( arma::vec y, arma::mat Sigma ) {

  int p = Sigma.n_rows;


  // inverse Sigma
  arma::mat Sigma1 = arma::inv_sympd(Sigma);
  // determinant Sigma
  double det_Sigma = arma::log_det_sympd(Sigma);
  // distance
  double dist = dmvnorm_distance( y, Sigma1);
  double pi1 = 3.14159265358979;
  double l1 = - p * std::log(2*pi1) - dist - det_Sigma;
  double ll = 0.5 * l1;
  return ll;
}
