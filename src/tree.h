#ifndef GUARD_tree_h
#define GUARD_tree_h

#include <map>
#include <cmath>
#include <cstddef>
#include "common.h"

//--------------------------------------------------
//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
//left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

//--------------------------------------------------
//info contained in a node, used by input operator
struct node_info {
  std::size_t id; //node id
  std::size_t v;  //variable
  std::size_t c;  //cut point
  double theta;   //theta
};

class tree {
public:
  //friends--------------------
  friend std::istream& operator>>(std::istream&, tree&);
  typedef tree* tree_p;
  typedef const tree* tree_cp;
  typedef std::vector<tree_p> npv;
  typedef std::vector<tree_cp> cnpv;
  tree();
  tree(const tree& n): theta(0.0),v(0),c(0),p(0),l(0),r(0) {cp(this,&n);}
  tree(double itheta): theta(itheta),v(0),c(0),p(0),l(0),r(0) {}
  void tonull(); //like a "clear", null tree has just one node
  ~tree() {tonull();}
  //operators----------
  tree& operator=(const tree&);
  tree& operator<=(tree&);
  //interface----------
  //set
  void settheta(double theta) {this->theta=theta;}
  void setv(size_t v) {this->v = v;}
  void setc(size_t c) {this->c = c;}
  //get
  double gettheta() const {return theta;}
  size_t getv() const {return v;}
  size_t getc() const {return c;}
  tree_p getp() {return p;}
  tree_p getl() {return l;}
  tree_p getr() {return r;}
  //tree functions--------------------
  tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there
  void pr(bool pc=true); //to screen, pc is "print children"
  //void pr(bool pc=true); //to screen, pc is "print children"
  size_t treesize(); //number of nodes in tree
  size_t nnogs();    //number of nog nodes (no grandchildren nodes)
  size_t nbots();    //number of bottom nodes
  bool birth(size_t nid, size_t v, size_t c, double thetal, double thetar);
  bool death(size_t nid, double theta);
  void birthp(tree_p np,size_t v, size_t c, double thetal, double thetar);
  void deathp(tree_p nb, double theta);
  // void updatesingletree(Rcpp::DataFrame df);
  void getbots(npv& bv);         //get bottom nodes
  void getnogs(npv& nv);         //get nog nodes (no granchildren)
  void getnodes(npv& v);         //get vector of all nodes
  void getnodes(cnpv& v) const;  //get vector of all nodes (const)
  tree_p bn(double *x,xinfo& xi); //find Bottom Node
  //std::vector<size_t> getbotsid(double *x, xinfo& xi, std::vector<size_t>& idv); // get bottom node id for a dataset
  void rg(size_t v, int* L, int* U); //recursively find region [L,U] for var v
  //node functions--------------------
  size_t nid() const; //nid of a node
  size_t depth();  //depth of a node
  char ntype(); //node type t:top, b:bot, n:no grandchildren i:interior (t can be b)
  bool isnog();
  void printtree(xinfo& xi);
  void printtreeToFile(xinfo& xi);
  size_t getbadcut(size_t v);
#ifndef NoRcpp
  Rcpp::List tree2list(xinfo& xi, double center=0., double scale=1.); // create an efficient list from a single tree
  //tree list2tree(Rcpp::List&, xinfo& xi); // create a tree from a list and an xinfo
  Rcpp::IntegerVector tree2count(size_t nvar); // for one tree, count the number of branches for each variable
  void initialize(size_t m);
#endif
private:
  double theta; //univariate double parameter
  //rule: left if x[v] < xinfo[v][c]
  size_t v;
  size_t c;
  //tree structure
  tree_p p; //parent
  tree_p l; //left child
  tree_p r; //right child
  //utility functions
  void cp(tree_p n,  tree_cp o); //copy tree
};
std::istream& operator>>(std::istream&, tree&);
std::ostream& operator<<(std::ostream&, const tree&);
#endif
