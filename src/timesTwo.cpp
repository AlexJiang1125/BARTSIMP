//#include <Rcpp.h>
#include "tree.h"
#include "rn.h"
#include <iostream>

//using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' Multiply a number by two
//'
//' @useDynLib BARTSIMP, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesTwo() {
  tree mytree;
  std::vector<double> v = {1.0,2.0,3.0,4.0,5.0};
  double r = log_sum_exp(v);
  cout << "logsumexp is " << r << endl;
  cout << "tree size is " << mytree.treesize() << endl;
  cout << "node is " << mytree.getp() << endl;
  mytree.pr(true);
  return 0;
}


