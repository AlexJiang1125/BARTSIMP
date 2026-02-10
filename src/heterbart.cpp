/*
 */

#define ARMA_USE_SUPERLU 1
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace std::chrono;
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
double heterbart::draw_sigmaupdate(double *sigma, rn& gen, double nu, double lambda, double &kappa, double &sigma_m, double &mlik, bool isexact)
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
    }
    // standard T
    heterbd_test(t[j],xi,r,di,pi,sigma,nv,pv,aug,gen,sigma_m,kappa,isexact);
    //hetergetmargprob_test(t[j], xi, di, r, sigma_m, pi);
    //hetergetmargprob_testonly(t[j], xi, di, r, sigma_m, pi);
    // standard mu
    //cout << "spatvars " << sigma_m << " " << kappa << " nonspatvars " << sigma[0] << endl;
    //auto t1 = high_resolution_clock::now();

    heterdrmu(t[j],xi,di,pi,sigma,gen);

    //heterdrmu_new (t[j],xi,di,pi,sigma[0],gen,r,sigma_m, kappa);
    //auto t2 = high_resolution_clock::now();
    //auto ms_int = duration_cast<milliseconds>(t2 - t1);
    //cout << "origin: " << ms_int.count() << endl;
    //std::ofstream MyFile("/Users/alexziyujiang/Documents/GitHub/BART_SIMP/results/results_origin_largedata.txt", MyFile.out | MyFile.app);
    //MyFile << ms_int.count() << " ";
    //MyFile << endl;
    //auto t11 = high_resolution_clock::now();
    //heterdrmu_new_approx(t[j],xi,di,pi,sigma[0],gen,r,sigma_m, kappa);
    //auto t22 = high_resolution_clock::now();
    //auto ms_int2 = duration_cast<milliseconds>(t22 - t11);
    //cout << "approx: " << ms_int2.count() << endl;
    //std::ofstream MyFile2("/Users/alexziyujiang/Documents/GitHub/BART_SIMP/results/results_approx_largedata.txt", MyFile2.out | MyFile2.app);
    //MyFile2 << ms_int2.count() << " ";
    //MyFile2 << endl;

    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += ftemp[k];
      //cout << "(" << ftemp[k] << "," << allfit[k] << ") ";
    }
  }
  double r_sigma[n];
  Rcpp::NumericVector r_sigma_(n);
  for(size_t k=0;k<n;k++) {
    r_sigma[k] = y[k] - allfit[k];
  }
  for(size_t k=0; k<n; k++) {
    r_sigma_[k] = r_sigma[k];
  }
  //cout << "sd " << Rcpp::sd(r_sigma_);
  sigma[0] = heterbd_drawsigma(di, r_sigma, sigma[0],kappa, sigma_m, pi, nu, lambda, gen, isexact);
  heterbd_drawspathyperpars(di, r_sigma,sigma[0], kappa, sigma_m, pi, mlik, gen, isexact);

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
void heterbart::draw_test(double *sigma, rn& gen, double &kappa, double &sigma_m, bool isexact)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
      //cout << allfit[k] << " ";
    }
    heterbd_test(t[j],xi,r,di,pi,sigma,nv,pv,aug,gen,sigma_m,kappa,isexact);
    heterdrmu(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += ftemp[k];
    }
  }
}

// covert to sparse matrix
arma::sp_mat heterbart::convertSparse(Rcpp::S4 mat) {

  // obtain dim, i, p. x from S4 object
  IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

  int nrow = dims[0], ncol = dims[1];

  // use Armadillo sparse matrix constructor
  arma::sp_mat res(i, p, x, nrow, ncol);
  return res;
  // Rcout << "SpMat res:\n" << res << std::endl;
}


// create mesh
void heterbart::makemesh() {
  Environment env("package:INLA");
  Rcpp::List myList(2);
  Function inlaMesh2D = env["inla.mesh.2d"];
  // Rcpp::NumericVector edge_arg(2);
  // edge_arg[0] = 0.04;
  // edge_arg[1] = 0.08;
  //
  // Rcpp::NumericVector offset(2);
  // offset[0] = 0.04;
  // offset[1] = 0.08;

  // Rcpp::NumericVector domainvec(8);
  // domainvec[0] = 0.0;
  // domainvec[1] = 1.0;
  // domainvec[2] = 1.0;
  // domainvec[3] = 0.0;
  // domainvec[4] = 0.0;
  // domainvec[5] = 0.0;
  // domainvec[6] = 1.0;
  // domainvec[7] = 1.0;

  Rcpp::NumericVector edge_arg = {0.04, 0.08};//{0.04,0.08};
  Rcpp::NumericVector offset = {0.04,0.08};
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
  double edge_1 = 0.1*avg_lth;
  double edge_2 = 0.2*avg_lth;

  Rcpp::NumericVector edge_arg = {edge_1,edge_2};
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
  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  Rcpp::List myList_coords(2);
  Rcpp::CharacterVector namevec_coords;
  std::string namestem_coords = "Column Heading";
  myList_coords[0] = di.s1;
  namevec_coords.push_back("cx");
  myList_coords[1] = di.s2;
  namevec_coords.push_back("cy");
  myList.attr("names") = namevec_coords;
  Rcpp::DataFrame dfout_coords(myList_coords);
  NumericMatrix coords = internal::convert_using_rfunction(dfout_coords, "as.matrix");
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["loc"] = coords
  );
  di.A = Rcpp::as<arma::sp_mat>(A);
}

// create mesh, by points
void heterbart::makemesh_bypoints_unique() {
  Environment env("package:INLA");
  Rcpp::List myList(2);
  Function inlaMesh2D = env["inla.mesh.2d"];
  Rcpp::NumericVector domainvec = {0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0};
  domainvec.attr("dim") = Rcpp::Dimension(4,2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  Function unique("unique");
  Rcpp::NumericVector it_s1 = unique(di.s1);
  Rcpp::NumericVector it_s2 = unique(di.s2);

  myList[0] = it_s1;
  myList[1] = it_s2;
  double xmin = Rcpp::min(di.s1);
  double xmax = Rcpp::max(di.s1);
  double ymin = Rcpp::min(di.s2);
  double ymax = Rcpp::max(di.s2);
  double avg_lth = (ymax - ymin + xmax - xmin)*0.5;
  double edge_1 = 0.1*avg_lth;
  double edge_2 = 0.2*avg_lth;

  Rcpp::NumericVector edge_arg = {edge_1,edge_2};
  double cutoff = 0.05*avg_lth;
  namevec.push_back("cx");
  namevec.push_back("cy");
  myList.attr("names") = namevec;
  Rcpp::DataFrame dfout(myList);
  Rcpp::List mesh = inlaMesh2D(Rcpp::_["loc.domain"] = dfout,
                               Rcpp::_["max.edge"] = edge_arg,
                               Rcpp::_["cutoff"] = cutoff);
  di.mesh = mesh;
  int nmesh = mesh["n"];
  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  Rcpp::List myList_coords(2);
  Rcpp::CharacterVector namevec_coords;
  std::string namestem_coords = "Column Heading";
  myList_coords[0] = it_s1;
  namevec_coords.push_back("cx");
  myList_coords[1] = it_s2;
  namevec_coords.push_back("cy");
  myList.attr("names") = namevec_coords;
  Rcpp::DataFrame dfout_coords(myList_coords);
  NumericMatrix coords = internal::convert_using_rfunction(dfout_coords, "as.matrix");
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = di.mesh,
    Rcpp::_["loc"] = coords
  );
  di.A_unique = Rcpp::as<arma::sp_mat>(A);
}
