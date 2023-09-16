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

#include "heterbart.h"

//--------------------------------------------------
void heterbart::pr()
{
  cout << "+++++heterbart object:\n";
  bart::pr();
}
//--------------------------------------------------
Rcpp::List heterbart::printtrees()
{
  Rcpp::List res;
  size_t nn = t[0].treesize();
  for (size_t j=0; j<m; j++) {
    tree mytr = t[j];
    tree::cnpv nds;
    mytr.getnodes(nds);
    Rcpp::NumericVector nodeid;
    Rcpp::NumericVector splitvarid;
    Rcpp::NumericVector varcutid;
    Rcpp::NumericVector varcutval;
    Rcpp::NumericVector nodeval;
    for (size_t k=0; k<nds.size(); k++) {
      nodeid.push_back(nds[k]->nid());
      splitvarid.push_back(nds[k]->getv());
      varcutid.push_back(nds[k]->getc());
      varcutval.push_back(xi[nds[k]->getv()][nds[k]->getc()]);
      nodeval.push_back(nds[k]->gettheta());
    }
    Rcpp::DataFrame treedf = Rcpp::DataFrame::create(_["nodeid"] = nodeid,
                                                 _["splitvarid"] = splitvarid,
                                                 _["varcutid"] = varcutid,
                                                 _["varcutval"] = varcutval,
                                                 _["nodeval"] = nodeval);
    res.push_back(treedf);
  }
  return res;
}
//--------------------------------------------------
void heterbart::initializetrees() {
  for (size_t j=0;j<m;j++) {
    t[j].initialize(m);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k] + ftemp[k];
    }
  }
}


//--------------------------------------------------
// with sigma updated
double heterbart::draw_sigmaupdate(double *sigma, rn& gen, double nu, double lambda, double &kappa, double &sigma_m, double &mlik)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
    }
    // standard T
    heterbd_test(t[j],xi,r,di,pi,sigma,nv,pv,aug,gen,sigma_m,kappa);
    // standard mu
    //cout << "spatvars " << sigma_m << " " << kappa << " nonspatvars " << sigma[0] << endl;
    heterdrmu_new(t[j],xi,di,pi,sigma[0],gen,r,sigma_m, kappa);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += ftemp[k];
      //cout << "(" << ftemp[k] << "," << allfit[k] << ") ";
    }
  }
  double r_sigma[n];
  for(size_t k=0;k<n;k++) {
    r_sigma[k] = y[k] - allfit[k];
  }
  sigma[0] = heterbd_drawsigma(di, r_sigma, sigma[0],kappa, sigma_m, pi, nu, lambda);
  heterbd_drawspathyperpars(di, r_sigma,sigma[0], kappa, sigma_m, pi, mlik);

  //if(dartOn) {
  //  draw_s(nv,lpv,theta,gen);
  //  draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
  //  for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
  //}
  return sigma[0];
}

void heterbart::draw(double *sigma, rn& gen, double nu, double lambda)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
    }
    //heterbd_new_test(t[j],xi,r,di,pi,sigma[0],nv,pv,aug,gen);
    heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen,r);
    //cout << "tree ok" << endl;
    heterdrmu(t[j],xi,di,pi,sigma,gen);
    //cout << "node ok" << endl;
    //heterdrmu_new_test(t[j],xi,di,pi,sigma[0],gen,r);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += ftemp[k];
    }
  }
  //if(dartOn) {
  //  draw_s(nv,lpv,theta,gen);
  //  draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
  //  for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
  //}
  return;
}



//--------------------------------------------------
void heterbart::updatetrees(Rcpp::List treeprev) {
  for(size_t j=0; j<m; j++) {
    Rcpp::List df = treeprev[j];
    updatetrees_module(df, 1, int(j));
  }
}

void heterbart::updatetrees_module(Rcpp::List df, int rowid, int treeid) {
  Rcpp::IntegerVector nodeids = df["nodeid"];
  Rcpp::IntegerVector splitvarids = df["splitvarid"];
  Rcpp::IntegerVector varcutids = df["varcutid"];
  Rcpp::NumericVector nodevals = df["nodeval"];
  int leftchild = rowid*2;
  int index = findindex(nodeids, rowid);
  // find position in the vector
  int index_left = findindex(nodeids, leftchild);
  if (index_left == -1) {
    return;
  } else {
    int index_right = findindex(nodeids, leftchild+1);
    t[treeid].birth(nodeids[index], splitvarids[index], varcutids[index], nodevals[index_left], nodevals[index_right]);
    updatetrees_module(df, leftchild, treeid);
    updatetrees_module(df, leftchild + 1, treeid);
  }
}

int heterbart::findindex(Rcpp::IntegerVector vec, int val) {
  std::vector<int> _vec(vec.begin(), vec.end());
  // find position in the vector
  std::vector<int>::iterator it = std::find(_vec.begin(), _vec.end(), val);
  if (it != _vec.end()) {
    int index = std::distance(_vec.begin(), it);
    return index;
  } else {
    return -1;
  }
}



//--------------------------------------------------
void heterbart::draw_test(double *sigma, rn& gen, double &kappa, double &sigma_m)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
      //cout << allfit[k] << " ";
    }
    heterbd_test(t[j],xi,r,di,pi,sigma,nv,pv,aug,gen,sigma_m,kappa);
    heterdrmu(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += ftemp[k];
    }
  }
}

// create mesh
void heterbart::makemesh() {
  Environment env("package:INLA");
  Rcpp::List myList(2);
  Function inlaMesh2D = env["inla.mesh.2d"];
  Rcpp::NumericVector edge_arg = {0.1/2.5,0.2/2.5};
  Rcpp::NumericVector offset = {0.1/2.5,0.2/2.5};
  Rcpp::NumericVector domainvec = {0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0};
  domainvec.attr("dim") = Rcpp::Dimension(4,2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  Rcpp::List mesh = inlaMesh2D(Rcpp::_["loc.domain"] = domainvec,
                               //Rcpp::_["loc.domain"] = domainvec,
                               Rcpp::_["max.edge"] = edge_arg,
                               Rcpp::_["offset"] = offset);
  di.mesh = mesh;
  int nmesh = mesh["n"];
}

// create mesh, by points
void heterbart::makemesh_bypoints() {
  Environment env("package:INLA");
  Rcpp::List myList(2);
  Function inlaMesh2D = env["inla.mesh.2d"];
  Rcpp::NumericVector domainvec = {0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0};
  domainvec.attr("dim") = Rcpp::Dimension(4,2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = di.s1;
  double xmin = Rcpp::min(di.s1);
  double xmax = Rcpp::max(di.s1);
  double ymin = Rcpp::min(di.s2);
  double ymax = Rcpp::max(di.s2);
  double avg_lth = (ymax - ymin + xmax - xmin)*0.5;
  Rcpp::NumericVector edge_arg = {0.1*avg_lth,0.2*avg_lth};
  double cutoff = 0.05*avg_lth;
  namevec.push_back("cx");
  myList[1] = di.s2;
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  Rcpp::List mesh = inlaMesh2D(Rcpp::_["loc.domain"] = dfout,
                               Rcpp::_["max.edge"] = edge_arg,
                               Rcpp::_["cutoff"] = cutoff);
  di.mesh = mesh;
  int nmesh = mesh["n"];
}
