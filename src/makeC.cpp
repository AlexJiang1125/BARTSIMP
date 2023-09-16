//#include <Rcpp.h>
//using namespace Rcpp;


#include "makeC.h"
#include "treefuns.h"

//' @export
// [[Rcpp::export]]
SEXP makeC(std::vector<size_t> idv, size_t n) {
  // get a vector with unique elements in idv
  std::set<int> s;
  unsigned size = idv.size();
  for( unsigned i = 0; i < size; ++i ) s.insert( idv[i] );
  std::vector<size_t> idv_unique;
  idv_unique.assign(s.begin(),s.end());
  if (idv_unique.size() > 1) {
    Rcpp::NumericVector v (n,0);
    Rcpp::DataFrame df = DataFrame::create(Named(std::to_string(idv_unique[0])) = v);
    for (int i = 1; i < idv_unique.size(); i++) {
      Rcpp::NumericVector v1 (n,0);
      df.push_back(clone(v1), std::to_string(idv_unique[i]));
    }
    for (int i = 0; i < idv.size(); i++) {
      for (int j = 0; j < idv_unique.size(); j++) {
        if (idv[i] == idv_unique[j]) {
          Rcpp::NumericVector vtemp = df[j];
          vtemp[i] = 1;
        }
      }
    }
    return df;
  } else {
    Rcpp::NumericVector v (n,1);
    return v;
  }
}




void printC(Rcpp::DataFrame df) {
  Rcpp::NumericVector x0 = df[0];
  std::ofstream MyFile("/Users/alexziyujiang/Documents/data/SBART/results.txt", MyFile.out | MyFile.app);
  int n = x0.size();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < df.size(); j++) {
      Rcpp::NumericVector temp = df[j];
      MyFile << temp[i] << " ";
    }
    MyFile << endl;
  }
  MyFile.close();
}


//' @export
// [[Rcpp::export]]
arma::mat convertC(SEXP mat) {
  return as<arma::mat>(mat);
}

//[[Rcpp::export]]
NumericMatrix testDFtoNM(DataFrame x) {
  NumericMatrix y = internal::convert_using_rfunction(x, "as.matrix");
  return y;
}

//' @export
// [[Rcpp::export]]
void callPrint(RObject x) {
  Rcpp::print(x);             // will work on any SEXP object
}

// [[Rcpp::export]]
void useOperatorOnVector(NumericVector x) {
  Rcpp::Rcout << "" << std::endl << x << std::endl;
}

// [[Rcpp::export]]
void useOperatorOnMatrix(NumericMatrix x) {
  Rcpp::Rcout << "" << std::endl << x << std::endl;
}
