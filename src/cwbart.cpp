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

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cwbart(
    SEXP _in,            //number of observations in training data
    SEXP _ip,		//dimension of x
    SEXP _inunique, // number of unique observations in training data
    SEXP _inp,		//number of observations in test data
    SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _ixunique,
    SEXP _iy,		//y, train,  nx1
    SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _s1,
    SEXP _s2,
    SEXP _weighted,
    SEXP _size_dat,
    SEXP _im,		//number of trees
    SEXP _inc,		//number of cut points
    SEXP _ind,		//number of kept draws (except for thinnning ..)
    SEXP _iburn,		//number of burn-in draws skipped
    SEXP _ipower,
    SEXP _ibase,
    SEXP _itau,
    SEXP _inu,
    SEXP _ilambda,
    SEXP _isigest,
    SEXP _irho0,
    SEXP _isigmam0,
    SEXP _irange_new,
    SEXP _isig_m_new,
    SEXP _isigest_new,
    SEXP _iw,
    SEXP _idart,
    SEXP _itheta,
    SEXP _iomega,
    SEXP _igrp,
    SEXP _ia,
    SEXP _ib,
    SEXP _irho,
    SEXP _iaug,
    SEXP _inkeeptrain,
    SEXP _inkeeptest,
    SEXP _inkeeptestme,
    SEXP _inkeeptreedraws,
    SEXP _inprintevery,
    //   SEXP _treesaslists,
    SEXP _Xinfo,
    SEXP _iftest,
    SEXP _inwarmup,
    SEXP _usecoords,
    SEXP _treeupdate,
    SEXP _treeprev
)
{

  //--------------------------------------------------
  //process args
  size_t n = Rcpp::as<int>(_in);
  size_t p = Rcpp::as<int>(_ip);
  size_t np = Rcpp::as<int>(_inp);
  size_t nunique = Rcpp::as<int>(_inunique);
  Rcpp::NumericVector  xv(_ix);
  double *ix = &xv[0];
  Rcpp::NumericVector xunique(_ixunique);
  double *ixunique = &xunique[0];
  Rcpp::NumericVector  yv(_iy);
  double *iy = &yv[0];
  Rcpp::NumericVector sv1(_s1);
  double *s1 = &sv1[0];
  Rcpp::NumericVector sv2(_s2);
  double *s2 = &sv2[0];
  bool weighted = Rcpp::as<bool>(_weighted);
  Rcpp::NumericVector vsize_dat(_size_dat);
  double *size_dat = &vsize_dat[0];
  Rcpp::NumericVector  xpv(_ixp);
  double *ixp = &xpv[0];
  size_t m = Rcpp::as<int>(_im);
  Rcpp::IntegerVector _nc(_inc);
  int *numcut = &_nc[0];
  //size_t nc = Rcpp::as<int>(_inc);
  size_t nd = Rcpp::as<int>(_ind);
  size_t burn = Rcpp::as<int>(_iburn);
  double mybeta = Rcpp::as<double>(_ipower);
  double alpha = Rcpp::as<double>(_ibase);
  double tau = Rcpp::as<double>(_itau);
  double nu = Rcpp::as<double>(_inu);
  double rho_0 = Rcpp::as<double>(_irho0);
  double sigmam_0 = Rcpp::as<double>(_isigmam0);
  double range_new = Rcpp::as<double>(_irange_new);
  double sig_m_new = Rcpp::as<double>(_isig_m_new);
  double sigest_new = Rcpp::as<double>(_isigest_new);
  double lambda = Rcpp::as<double>(_ilambda);
  double sigma=Rcpp::as<double>(_isigest);
  Rcpp::NumericVector  wv(_iw);
  double *iw = &wv[0];
  bool dart;
  if(Rcpp::as<int>(_idart)==1) dart=true;
  else dart=false;
  double a = Rcpp::as<double>(_ia);
  double b = Rcpp::as<double>(_ib);
  double rho = Rcpp::as<double>(_irho);
  bool aug;
  bool iftest = Rcpp::as<bool>(_iftest);
  bool usecoords = Rcpp::as<bool>(_usecoords);
  if(Rcpp::as<int>(_iaug)==1) aug=true;
  else aug=false;
  double theta = Rcpp::as<double>(_itheta);
  double omega = Rcpp::as<double>(_iomega);
  Rcpp::IntegerVector _grp(_igrp);
  int *grp = &_grp[0];
  size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
  size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
  size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
  size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
  size_t printevery = Rcpp::as<int>(_inprintevery);
  size_t nwarmup = Rcpp::as<int>(_inwarmup);
  //   int treesaslists = Rcpp::as<int>(_treesaslists);
  Rcpp::NumericMatrix Xinfo(_Xinfo);
  //   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

  //return data structures (using Rcpp)
  Rcpp::NumericVector trmean(n); //train
  Rcpp::NumericVector temean(np);
  Rcpp::NumericVector sdraw(nd+burn+nwarmup);
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);
  //   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
  Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
  Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
  Rcpp::IntegerMatrix lfcounts(nkeeptreedraws,p);
  // storing spatial hyperpars
  Rcpp::NumericVector kappas(nd+burn+nwarmup);
  Rcpp::NumericVector sigmams(nd+burn+nwarmup);
  // for debug use: log likelihood
  Rcpp::NumericVector mliks(nd+burn+nwarmup);
  Rcpp::List treeprev(_treeprev);
  bool treeupdate = Rcpp::as<bool>(_treeupdate);



  //random number generation
  arn gen;

  heterbart bm(m);

  // setinfo
  if(Xinfo.size()>0) {
    xinfo _xi;
    _xi.resize(p);
    for(size_t i=0;i<p;i++) {
      _xi[i].resize(numcut[i]);
      //Rcpp::IntegerVector cutpts(Xinfo[i]);
      for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
    }
    bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

  void cwbart(
      size_t n,            //number of observations in training data
      size_t p,		//dimension of x
      size_t nunique,
      size_t np,		//number of observations in test data
      double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
      double* iy,		//y, train,  nx1
      double* ixunique,
      double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
      double* s1,   //s1, nx1
      double* s2,   //s2, nx1,
      double* size_dat, // observation per cluster (used for scaling)
      size_t m,		//number of trees
      int* numcut,		//number of cut points
      size_t nd,		//number of kept draws (except for thinnning ..)
      size_t burn,		//number of burn-in draws skipped
      double mybeta,
      double alpha,
      double tau,
      double nu,
      double lambda,
      double sigma,
      double rho_0,
      double sigmam_0,
      double range_new,
      double sig_m_new,
      double sigest_new,
      double* iw,
      bool dart,
      double theta,
      double omega,
      int *grp,
      double a,
      double b,
      double rho,
      bool aug,
      size_t nkeeptrain,
      size_t nkeeptest,
      size_t nkeeptestme,
      size_t nkeeptreedraws,
      size_t printevery,
      //   int treesaslists,
      unsigned int n1, // additional parameters needed to call from C++
      unsigned int n2,
      double* trmean,
      double* temean,
      double* sdraw,
      double* _trdraw,
      double* _tedraw,
      bool iftest,
      bool usecoords,
      size_t nwarmup,
      bool treeupdate,
      Rcpp::List treeprev
  )
  {

    //return data structures (using C++)
    std::vector<double*> trdraw(nkeeptrain);
    std::vector<double*> tedraw(nkeeptest);

    for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
    for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

    std::vector< std::vector<size_t> > varcnt;
    std::vector< std::vector<double> > varprb;

    //random number generation
    arn gen(n1, n2);

    heterbart bm(m);

#endif

    for(size_t i=0;i<n;i++) trmean[i]=0.0;
    for(size_t i=0;i<np;i++) temean[i]=0.0;

    printf("*****Into main of wbart\n");
    //-----------------------------------------------------------

    size_t skiptr,skipte,skipteme,skiptreedraws;
    if(nkeeptrain) {skiptr=nd/nkeeptrain;}
    else skiptr = nd+1;
    if(nkeeptest) {skipte=nd/nkeeptest;}
    else skipte=nd+1;
    if(nkeeptestme) {skipteme=nd/nkeeptestme;}
    else skipteme=nd+1;
    if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
    else skiptreedraws=nd+1;

    //--------------------------------------------------
    //print args
    lambda = 5;
    nu = 5;

    printf("*****Data:\n");
    printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
    printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
    printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
    if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
    printf("*****Number of Trees: %zu\n",m);
    printf("*****Size: %lf, %lf\n", size_dat[0],size_dat[nunique-1]);
    printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
    printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
    printf("*****Prior:beta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
           mybeta,alpha,tau,nu,lambda);
    printf("*****Pcprior:rho0, sigmam0: %lf, %lf\n", rho_0, sigmam_0);
    printf("*****Spathyper init points :rho, sigmam, sigma: %lf, %lf, %lf\n", range_new, sig_m_new, sigest_new);
    printf("*****sigma: %lf\n",sigma);
    printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
    cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: "
         << dart << ',' << theta << ',' << omega << ',' << a << ','
         << b << ',' << rho << ',' << aug << endl;
    printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
           nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
    printf("*****printevery: %zu\n",printevery);
    printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);

    //--------------------------------------------------
    //heterbart bm(m);
    bm.setprior(alpha,mybeta,tau);
    Rcpp::NumericVector rs1, rs2, rsize_dat;
    for (int i = 0; i < n; i++) {
      rs1.push_back(s1[i]);
      rs2.push_back(s2[i]);
    }
    for (int i = 0; i < nunique; i++) {
      rsize_dat.push_back(size_dat[i]);
    }
    bm.setdata(p,n,nunique,ix,iy,ixunique,rs1,rs2,rsize_dat,numcut[0]);
    bm.setdart(a,b,rho,aug,dart,theta,omega);
    cout << "usecoords" << usecoords;
    if (usecoords) {
      cout << "usecoords" << endl;
      bm.makemesh_bypoints();
    } else {
      cout << "notusecoords" << endl;
      bm.makemesh();
    }

    //--------------------------------------------------
    //sigma
    //gen.set_df(n+nu);
    double *svec = new double[n];
    for(size_t i=0;i<n;i++) svec[i]=iw[i]*sigma;

    // matern priors
    double mlik = -1000;
    bm.getpinfo().rho_0 = rho_0;
    bm.getpinfo().sigma_m0 = sigmam_0;
    bm.getpinfo().sigma2e = sigest_new*sigest_new;
    bm.getpinfo().sigma2x = sig_m_new*sig_m_new;
    bm.getpinfo().kappa = sqrt(8)/range_new;
    bm.getpinfo().alpha_1 = 0.05;
    bm.getpinfo().alpha_2 = 0.05;
    double kappa = bm.getpinfo().kappa;
    double sigma_m = std::sqrt(bm.getpinfo().sigma2x);


    //double& sigma_m = std::sqrt(sigma_m2);
    //double& sigma_m = std::pow(bm.getpinfo().sigma2x, 0.5);

    //--------------------------------------------------

    std::stringstream treess;  //string stream to write trees to
    treess.precision(10);
    treess << nkeeptreedraws << " " << m << " " << p << endl;
    // dart iterations
    std::vector<double> ivarprb (p,0.);
    std::vector<size_t> ivarcnt (p,0);

    //--------------------------------------------------
    //temporary storage
    //out of sample fit
    double* fhattest=0; //posterior mean for prediction
    if(np) { fhattest = new double[np]; }
    double restemp=0.0,rss=0.0;


    //--------------------------------------------------
    //mcmc
    printf("\nMCMC\n");
    //size_t index;
    size_t trcnt=0; //count kept train draws
    size_t tecnt=0; //count kept test draws
    size_t temecnt=0; //count test draws into posterior mean
    size_t treedrawscnt=0; //count kept bart draws
    bool keeptest,keeptestme,keeptreedraw;

    time_t tp;
    int time1 = time(&tp);
    xinfo& xi = bm.getxinfo();
    Rcpp::List treeout;
    bm.initializetrees();
    cout << "treeupdate " << treeupdate << endl;
    if (treeupdate) {
      bm.updatetrees(treeprev);
    }
    // std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt");
    // MyFile << "start" << endl;
    // MyFile.close();

    //cout << (iftest) << endl;

    for(size_t i=0;i<(nd+nwarmup+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      //if(i< nwarmup) {
      //  iftest = false;
      //} else {
      //  iftest = true;
      //}
      if (i == nwarmup) {
        cout << "warm up period ends!" << endl;
      }
      if (iftest == true) {
        cout << "sigma " << svec[0] << endl;
        bm.draw_sigmaupdate(svec,gen, nu, lambda, kappa, sigma_m, mlik);
        mliks[i] = mlik;
        kappas[i] = kappa;
        sigmams[i] = sigma_m;
        //bm.getpinfo().kappa = kappa;
        //bm.getpinfo().sigma2x = sigma_m*sigma_m;
        //bm.getpinfo().sigma2e = svec[0]*svec[0];
        //cout << kappa<<  ' ' << sigma_m << endl;
      } else {
        bm.draw(svec,gen, nu, lambda);
      }
      if(i % 1 ==0) {
        // output tree structure
        treeout.push_back(bm.printtrees());
      }
      //draw sigma
      if (iftest == false) {
        double flag = 0; // if flag 1, MH, else gibbs
        if (flag == 1) {
          for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k))/(iw[k]); rss += restemp*restemp;}
          double sigma_0 = svec[0];
          //double rands = gen.uniform();
          double rands = arma::randn<double>();
          double temp = rands*(0.2)+(2.0)*std::log(sigma_0);
          //double sigma_star = rands*(0.1)+sigma_0 - 0.05;
          double sigma_star = std::sqrt(std::exp(temp));
          //cout << sigma_star << endl;

          if (sigma_star < 0) {
            svec[0] = sigma;
            sdraw[i] = svec[0];
          } else {

          double lik1 = -rss/2.0/sigma_star/sigma_star-n/2.0*std::log(sigma_star)*(2.0);
          double lik2 = -rss/2.0/sigma_0/sigma_0-n/2.0*std::log(sigma_0)*(2.0);
          bool flag_true = 1;
          double p1 = R::dchisq(nu*lambda/sigma_star/sigma_star,nu,flag_true);
          double p2 = R::dchisq(nu*lambda/sigma/sigma,nu,flag_true);
          double mh_ratio = lik1 - lik2 + (p1 - p2) + (std::log(sigma_0) - std::log(sigma_star))*(4.0) - (std::log(sigma_0) - std::log(sigma_star))*(2.0);
          if (mh_ratio > 0) {
            svec[0] = sigma_star;
            sdraw[i] = svec[0];
          } else {
            double random_gen = gen.uniform();
            if (random_gen < std::exp(mh_ratio)) {
              svec[0] = sigma_star;
              sdraw[i] = svec[0];
            }
            svec[0] = sigma_0;
            sdraw[i] = svec[0];
          }
          }
        } else {

        // OG parts
        for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k))/(iw[k]); rss += restemp*restemp;}
        sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
        for(size_t k=0;k<n;k++) svec[k]=iw[k]*sigma;
        sdraw[i]=sigma;
        mliks[i] = -100000;
        kappas[i] = kappa;
        sigmams[i] = sigma_m;
        }
        // end OG
      } else {
        sdraw[i]=svec[0];
        //cout << iy[0]-bm.f(0) << "true residual" << endl;
      }
      rss=0.0;
      cout << "startprint" << endl;
      if(i>=burn+nwarmup) {
        //size_t* bnds = bm.getLeafNodes();
        //for (size_t k=0;k<p;k++) {
        //  lfcounts[i-burn+1,k] = bnds[i-burn];
        //}
        // sum of tree fitted values
        // to be divided over iterations later
        for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);

        if(nkeeptrain && (((i-(burn+nwarmup)+1) % skiptr) ==0)) {
          //index = trcnt*n;;
          //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
          cout << "trcnt" << trcnt << endl;
          for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
          trcnt+=1;
        }
        keeptest = nkeeptest && (((i-(burn+nwarmup)+1) % skipte) ==0) && np;
        keeptestme = nkeeptestme && (((i-(burn+nwarmup)+1) % skipteme) ==0) && np;
        if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
        if(keeptest) {
          //index=tecnt*np;
          //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
          for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
          tecnt+=1;
        }
        if(keeptestme) {
          for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
          temecnt+=1;
        }
        keeptreedraw = nkeeptreedraws && (((i-(burn+nwarmup)+1) % skiptreedraws) ==0);
        if(keeptreedraw) {
          //	   #ifndef NoRcpp
          //	   Rcpp::List lists(m*treesaslists);
          //	   #endif

          for(size_t j=0;j<m;j++) {
            treess << bm.gettree(j);
            /*
#ifndef NoRcpp
             varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
             if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
#endif
             */
          }
#ifndef NoRcpp
          //	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
          ivarcnt=bm.getnv();
          ivarprb=bm.getpv();
          size_t k=(i-(burn+nwarmup))/skiptreedraws;
          //cout << "k: " << k << endl;
          for(size_t j=0;j<p;j++){
            varcnt(k,j)=ivarcnt[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb(k,j)=ivarprb[j];
            //varprb(i-burn,j)=ivarprb[j];
          }
#else
          varcnt.push_back(bm.getnv());
          varprb.push_back(bm.getpv());
#endif

          treedrawscnt +=1;
        }
      }
      cout << "endprint" << endl;
    }
    int time2 = time(&tp);
    printf("time: %ds\n",time2-time1);
    for(size_t k=0;k<n;k++) trmean[k]/=nd;
    for(size_t k=0;k<np;k++) temean[k]/=temecnt;
    printf("check counts\n");
    printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
    //--------------------------------------------------
    //PutRNGstate();

    if(fhattest) delete[] fhattest;
    if(svec) delete [] svec;

    //--------------------------------------------------
    //return
#ifndef NoRcpp
    Rcpp::List ret;
    ret["sigma"]=sdraw;
    ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    ret["yhat.test.mean"]=temean;
    ret["yhat.test"]=tedraw;
    //ret["varcount"]=varcount;
    ret["varcount"]=varcnt;
    ret["varprob"]=varprb;
    ret["treeout"]=treeout;
    ret["leafcounts"]=lfcounts;

    //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
    //}

    Rcpp::List xiret(xi.size());
    for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }

    Rcpp::List treesL;
    //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
    //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
    //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
    treesL["cutpoints"] = xiret;
    treesL["trees"]=Rcpp::CharacterVector(treess.str());
    //   if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws"] = treesL;
    // for (int i = 0; i < sigmams.size(); i++) {
    //   cout << sigmams[i] << " " << kappas[i] << " " << mliks[i] << endl;
    // }
    ret["sigmams"] = sigmams;
    ret["kappas"] = kappas;
    //ret["mliks"] = mliks;
    return ret;
#else

#endif

  }
