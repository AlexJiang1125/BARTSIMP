/*
 */


#include "bart.h"
using namespace arma;
//--------------------------------------------------
//constructor
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(const bart& ib):m(ib.m),t(m),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false)
{
  this->t = ib.t;
}
bart::~bart()
{
  if(allfit) delete[] allfit;
  if(r) delete[] r;
  if(ftemp) delete[] ftemp;
}

//--------------------------------------------------
//operators
bart& bart::operator=(const bart& rhs)
{
  if(&rhs != this) {

    this->t = rhs.t;
    this->m = t.size();

    this->pi = rhs.pi;

    p=0;n=0;x=0;y=0;
    xi.clear();

    if(allfit) {delete[] allfit; allfit=0;}
    if(r) {delete[] r; r=0;}
    if(ftemp) {delete[] ftemp; ftemp=0;}

  }
  return *this;
}
//--------------------------------------------------
//get,set
void bart::setm(size_t m)
{
  t.resize(m);
  this->m = t.size();

  if(allfit && (xi.size()==p)) predict(p,n,x,allfit);
}

//--------------------------------------------------
void bart::setxinfo(xinfo& _xi)
{
  size_t p=_xi.size();
  xi.resize(p);
  for(size_t i=0;i<p;i++) {
    size_t nc=_xi[i].size();
    xi[i].resize(nc);
    for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
  }
}

//--------------------------------------------------
void bart::setdata(size_t p, size_t n, size_t nunique, double *x, double *y,double *xunique,
                   Rcpp::NumericVector s1, Rcpp::NumericVector s2,
                   Rcpp::NumericVector dat_size, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, nunique, x, y, xunique, nc,s1, s2, dat_size);
  delete [] nc;
}

void bart::setdata(size_t p, size_t n, size_t nunique, double *x, double *y, double *xunique, int *nc,
                   Rcpp::NumericVector s1, Rcpp::NumericVector s2,
                   Rcpp::NumericVector dat_size)
{
  this->p=p; this->n=n; this->x=x; this->y=y;
  cout << "xi size is : " << xi.size() << endl;
  if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

  if(allfit) delete[] allfit;
  allfit = new double[n];
  predict(p,n,x,allfit);

  if(r) delete[] r;
  r = new double[n];

  if(ftemp) delete[] ftemp;
  ftemp = new double[n];

  di.n=n; di.p=p; di.x = &x[0]; di.y=r;di.nunique = nunique;
  di.xunique = &xunique[0];
  di.s1 = s1; di.s2 = s2; di.dat_size = dat_size;

  for(size_t j=0;j<p;j++){
    nv.push_back(0);
    pv.push_back(1/(double)p);
  }
}
//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
  //uses: m,t,xi
{
  double *fptemp = new double[n];

  for(size_t j=0;j<n;j++) fp[j]=0.0;
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,fptemp);
    for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
  }

  delete[] fptemp;
}
//--------------------------------------------------
void bart::draw(double sigma, rn& gen)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      r[k] = y[k]-allfit[k];
    }
    bd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
    drmu(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
  }
  if(dartOn) {
    draw_s(nv,lpv,theta,gen);
    draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
    for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
  }
}
//--------------------------------------------------
//public functions
void bart::pr() //print to screen
{
  cout << "*****bart object:\n";
  cout << "m: " << m << std::endl;
  cout << "t[0]:\n " << t[0] << std::endl;
  cout << "t[m-1]:\n " << t[m-1] << std::endl;
  cout << "prior and mcmc info:\n";
  pi.pr();
  if(dart){
    cout << "*****dart prior (On):\n";
    cout << "a: " << a << std::endl;
    cout << "b: " << b << std::endl;
    cout << "rho: " << rho << std::endl;
    cout << "augmentation: " << aug << std::endl;
  }
  else cout << "*****dart prior (Off):\n";
  if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
  else cout << "data not set\n";
}

size_t* bart::getLeafNodes() {
  size_t bnds[m];
  for (size_t i=0; i<m; i++) {
    bnds[i] = t[i].nbots();
  }
  return bnds;
}

