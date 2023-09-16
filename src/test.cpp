#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bd.h"
#include "bart.h"
#include "heterbartfuns.h"
#include "heterbd.h"
#include "heterbart.h"
#include "makeC.h"
#include "rn.h"




//' @export
// [[Rcpp::export]]
double test(Rcpp::NumericVector s1,
         Rcpp::NumericVector s2,
         Rcpp::NumericVector ydata) {
  /*
   N: size of training data
   M: size of testing data
   P: covariates
   H: trees
   D: MCMC draws
   B: burn-in discarded
   C: cutpoints (an array since some variables may be discrete)
   */
  // H: number of trees
  size_t N=100, M=1, P=2, H=5, B=0, D=1, T=B+D; // B=250, D=1000,
  int C=100, augment=0, sparse=0, theta=0, omega=1;
  unsigned int n1=111;
  unsigned int n2=222; // seeds
  double k=2., power=2., base=0.95, a=0.5, b=1., tau=14./(2.*k*sqrt(H)),
    nu=3., ymu=0., ysd=0., sigquant=0.9, lambda, sigma=1.;
  Rcpp::NumericVector _xtrain(N*P);
  Rcpp::NumericVector _y(N);
  _xtrain[Rcpp::seq(0,N*P-1)] = Rcpp::runif(N*P,0,1);
  _y[Rcpp::seq(0,N-1)] = Rcpp::runif(N,0,1);
  double xtrain[N*P];
  double xunique[N*P];
  double y[N];
  for (int i=0;i<N*P; i++) {
    xtrain[i] = _xtrain[i];
    xunique[i] = _xtrain[i];
  }
  for (int i=0;i<N; i++) {
    y[i] = _y[i];
  }
  //double xtrain[100]= Rcpp::runif(100,0,1);
  //double xtest[2]={0., 1.};
  //double y[10]={0., 1., -1., 2., -2., 10., 11., 9., 12., 8.}; // outcome
  //double w[10]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; // known weights for SDi

  // tune lambda parameter
  for(size_t i=0; i<N; ++i) ymu += y[i]/N;
  for(size_t i=0; i<N; ++i) ysd += pow(y[i]-ymu, 2.)/(N-1.);
  ysd=sqrt(ysd);

  double* sd=new double[B+D];
  double* _ftrain=new double[N*D];
  double* _ftest =new double[M*D];

  std::vector<double*> ftrain(D);
  std::vector<double*> ftest(D);

  for(size_t h=0; h<D; ++h) ftrain[h]=&_ftrain[h*N];
  for(size_t h=0; h<D; ++h) ftest[h] =&_ftest[h*M];

  // test sample five trees
  heterbart Alg(H);
  Alg.setprior(base, power, tau);
  Alg.setdata(P, N, N, xtrain, y,xunique, s1,s2,C);
  Alg.setdart(a, b, P, augment, sparse, theta, omega);
  // print information for the first tree
  tree& t1 = Alg.gettree(0);
  cout << t1 << endl;
  xinfo& xif1 = Alg.getxinfo();
  // grow a tree
  tree::tree_p nx; //bottom node
  size_t v,c; //variable and cutpoint
  double pr; //part of metropolis ratio from proposal and prior
  tree::npv goodbots;  //nodes we could birth at (split on)
  pinfo& pi = Alg.getpinfo();
  cout << "pinfo: " << pi.sigma2x << endl;
  double PBx = getpb(t1,xif1,pi,goodbots); //prob of a birth at x
  cout << "pbx: " << PBx << endl;
  arn gen; // abstract random number generator
  bool aug = false;
  tree::npv bots;
  t1.getbots(bots);
  // t1.birthp(bots[0],0,20,0.2,0.6);
  // cout << t1 << endl;
  // tree::npv bots2;
  // t1.getbots(bots2);
  // t1.birthp(bots2[0],1,15,0.1,0.4);
  // cout << t1 << endl;
  std::vector<size_t> idv1;
  getbotsid(t1, xtrain, xif1, N, P, idv1);
  for (int i=0; i<idv1.size(); i++) {
    cout << idv1[i];
  }
  cout << endl;
  std::set<int> s;
  unsigned size = idv1.size();
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv2;
  idv2.assign( s.begin(), s.end() );
  for (int i=0; i<idv2.size(); i++) {
    cout << idv2[i];
  }
  Rcpp::DataFrame df = makeC(idv1, N);
  for (int i = 0; i < idv2.size(); i++) {
    useOperatorOnVector(df[i]);
  }
  Rcpp::CharacterVector cnames = df.names();
  for (int i = 0; i < cnames.size(); i++) {
    cout << cnames[i];
  }
  cout << endl;
  Environment env("package:INLA");
  Function inla = env["inla"];
  Rcpp::List myList(2);
  Rcpp::CharacterVector namevec;
  std::string namestem = "Column Heading";
  myList[0] = s1;
  namevec.push_back("cx");
  myList[1] = s2;
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
  Function inlaMesh2D = env["inla.mesh.2d"];
  Rcpp::NumericVector edge_arg = {1.0,1.0};
  Rcpp::List mesh = inlaMesh2D(Rcpp::_["loc"] = dfout,
                               Rcpp::_["max.edge"] = edge_arg);
  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];
  Rcpp::NumericVector priorRange = {10.0,0.5};
  Rcpp::NumericVector priorSigma = {0.5,0.5};
  Rcpp::List spde = inlaSpde2PcMatern(
    Rcpp::_["mesh"] = mesh,
    Rcpp::_["prior.range"] = priorRange,
    Rcpp::_["prior.sigma"] = priorSigma
  );
  Function inlaSpdeMakeA = env["inla.spde.make.A"];
  NumericMatrix coords = internal::convert_using_rfunction(dfout, "as.matrix");
  Rcpp::S4 A = inlaSpdeMakeA(
    Rcpp::_["mesh"] = mesh,
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
    Named("A") = inlaStackA(stkDat),
    Named("compute") = Rcpp::LogicalVector(TRUE)
  );
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula1),
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = Rcpp::LogicalVector(TRUE)),
                                 Rcpp::_["data"] = inlaStackData(stkDat, Rcpp::_["spde"] = spde),
                                 Rcpp::_["control.predictor"] = controlpred_list);
  Rcpp::NumericVector vout = inla_results["mlik"]; //storing the vector of mliks
  //Rcout << vout[0] << std::endl;

  return vout[0];
}


//' @export
// [[Rcpp::export]]
double test2(Rcpp::DataFrame mydf) {
  /*
   N: size of training data
   M: size of testing data
   P: covariates
   H: trees
   D: MCMC draws
   B: burn-in discarded
   C: cutpoints (an array since some variables may be discrete)
   */
  // H: number of trees
  Rcpp::NumericVector ydata = mydf[0];
  Rcpp::NumericVector s1 = mydf[1];
  Rcpp::NumericVector s2 = mydf[2];
  size_t N = mydf.nrows();
  size_t P = mydf.size() - 3;
  cout << P << endl;
  // load in data
  double xtrain[N*P];
  double xunique[N*P];
  double y[N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < P; j++) {
      Rcpp::NumericVector temp = mydf[j+3];
      xtrain[i*P+j] = temp[i];
      xunique[i*P+j] = temp[i];
    }
    y[i] = ydata[i];
  }
  size_t M=1, H=5, B=0, D=1, T=B+D; // B=250, D=1000,
  int C=100, augment=0, sparse=0, theta=0, omega=1;
  unsigned int n1=111;
  unsigned int n2=222; // seeds
  double k=2., power=2., base=0.95, a=0.5, b=1., tau=14./(2.*k*sqrt(H)),
    nu=3., ymu=0., ysd=0., sigquant=0.9, lambda, sigma=0.05;
  // the xtrains are randomly generated
  //Rcpp::NumericVector _xtrain(N*P);
  //Rcpp::NumericVector _y(N);
  //_xtrain[Rcpp::seq(0,N*P-1)] = Rcpp::runif(N*P,0,1);
  //_y[Rcpp::seq(0,N-1)] = Rcpp::runif(N,0,1);
  //for (int i=0;i<N*P; i++) {
  //  xtrain[i] = _xtrain[i];
  //}
  //for (int i=0;i<N; i++) {
  //  y[i] = _y[i];
  //}
  //double xtrain[100]= Rcpp::runif(100,0,1);
  //double xtest[2]={0., 1.};
  //double y[10]={0., 1., -1., 2., -2., 10., 11., 9., 12., 8.}; // outcome
  //double w[10]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; // known weights for SDi

  // tune lambda parameter
  for(size_t i=0; i<N; ++i) ymu += y[i]/N;
  for(size_t i=0; i<N; ++i) ysd += pow(y[i]-ymu, 2.)/(N-1.);
  ysd=sqrt(ysd);

  double* sd=new double[B+D];
  double* _ftrain=new double[N*D];
  double* _ftest =new double[M*D];

  std::vector<double*> ftrain(D);
  std::vector<double*> ftest(D);

  for(size_t h=0; h<D; ++h) ftrain[h]=&_ftrain[h*N];
  for(size_t h=0; h<D; ++h) ftest[h] =&_ftest[h*M];

  // test sample five trees
  heterbart Alg(H);
  Alg.setprior(base, power, tau);
  Alg.setdata(P, N,N, xtrain, y, xunique, s1,s2,C);
  Alg.setdart(a, b, P, augment, sparse, theta, omega);
  // print information for the first tree
  tree& t1 = Alg.gettree(0);
  cout << t1 << endl;
  xinfo& xif1 = Alg.getxinfo();
  // grow a tree
  tree::tree_p nx; //bottom node
  size_t v,c; //variable and cutpoint
  double pr; //part of metropolis ratio from proposal and prior
  tree::npv goodbots;  //nodes we could birth at (split on)
  pinfo& pi = Alg.getpinfo();
  cout << "pinfo: " << pi.sigma2x << endl;

  double PBx = getpb(t1,xif1,pi,goodbots); //prob of a birth at x
  cout << "pbx: " << PBx << endl;
  arn gen; // abstract random number generator
  cout << gen.uniform() << " " << gen.uniform() <<" " << gen.uniform() << endl;
  bool aug = false;
  //tree::npv bots;
  // t1.getbots(bots);
  // t1.birthp(bots[0],0,50,0.2,0.6);
  // cout << t1 << endl;
  // tree::npv bots2;
  // t1.getbots(bots2);
  // t1.birthp(bots2[0],1,45,0.1,0.4);
  // cout << t1 << endl;
  // t1.printtree(xif1);
  // compute the fitted values
  double fitv[N];
  fit(t1, xif1, P, N, xtrain, fitv);
  double res[N];
  for (int i = 0; i < N; i++) {
    res[i] = y[i] - fitv[i];
  }
  dinfo& di = Alg.getdinfo();
  double *xx = di.x;
  std::vector<size_t> nv = Alg.getnv();
  std::vector<double> pv = Alg.getpv();
  heterbd_new(t1,xif1,res,di,pi,sigma,nv,pv,aug,gen);
  t1.printtree(xif1);
  std::vector<double> yhat_;
  for (int i = 0; i < di.n; i++) {
    yhat_.push_back(res[i]);
  }
  std::vector<size_t> idv1;
  arma::vec yhat(yhat_);
  getbotsid(t1, di.x, xif1, di.n, di.p, idv1);
  Rcpp::DataFrame df = makeC(idv1, di.n);
  Rcpp::NumericMatrix dfmat = testDFtoNM(df);
  arma::mat Cmat = convertC(dfmat);
  cout << "tau: " << pi.tau << " kappa: " << pi.kappa <<" nu: " << pi.nu << endl;
  //100 3
  //cout << Cmat.n_rows << " " << Cmat.n_cols << endl;
  arma::mat prec = calcPrecMat(di, pi, sigma);
  arma::mat temp1 = Cmat.t()*prec;
  arma::mat temp2(Cmat.n_cols , Cmat.n_cols, fill::eye);
  temp2 = temp2 * pi.tau;
  temp2 = (temp1*Cmat + temp2).i();
  arma::mat L = arma::chol(temp2, "lower");
  arma::vec meanvec = temp2*temp1*yhat;
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(123);
  arma::vec tempz(Cmat.n_cols, arma::fill::randn);
  arma::vec mudraws = meanvec + L*tempz;
  for (int i = 0; i < Cmat.n_cols; i++) {
    cout << mudraws[i] << " ";
  }
  cout << endl;
  std::set<int> s;
  unsigned size = idv1.size();
  for( unsigned i = 0; i < size; ++i ) s.insert( idv1[i] );
  std::vector<size_t> idv_unique;
  idv_unique.assign(s.begin(),s.end());
  for (int i = 0; i < idv_unique.size(); i++) {
    int btid = idv_unique[i];
    tree::tree_p np = t1.getptr(btid);
    np->settheta(mudraws[i]);
  }
  t1.printtree(xif1);
  fit(t1, xif1, P, N, xtrain, fitv);
  double kappa = 2;
  double sigmam = 1;
  heterdrmu_new(t1, xif1, di, pi, 0.05, gen, res, kappa, sigmam);
  double db = hetergetmargprob(
    t1, xif1, di, res, 0.05, pi
  );

  return 0;
}


