#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' A prototype function that runs the inla function and
//' returns the marginal log-likelihood
//' @return inla_results the marginal log likelihood
//' @useDynLib BARTSIMP, .registration = TRUE
//' @param s1 a vector containing the longitude data
//' @param s2 a vector containing the latitude data
//' @param y the response variable
//' @importFrom Rcpp sourceCpp
//' @return The marginal likelihood for the model
//' @export
// [[Rcpp::export]]
RcppExport SEXP mySPDE(Rcpp::NumericVector s1,
                                      Rcpp::NumericVector s2,
                                      Rcpp::NumericVector y){
  // Obtaining namespace of INLA package
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
  std::string formula = "y ~ -1 + f(i, model = spde)";
  Function form("as.formula");
  Function inlaMesh2D = env["inla.mesh.2d"];
  Rcpp::NumericVector edge_arg = {1.0,1.0};
  Rcpp::List mesh = inlaMesh2D(Rcpp::_["loc"] = dfout,
                               Rcpp::_["max.edge"] = edge_arg);
  Function inlaSpde2PcMatern = env["inla.spde2.pcmatern"];
  Rcpp::NumericVector priorRange = {10.0,NA_REAL};
  Rcpp::NumericVector priorSigma = {0.5,NA_REAL};
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
  Rcpp::List y_list = Rcpp::List::create(Named("y") = y);
  Rcpp::List A_list = Rcpp::List::create(
    Named("A") = A,
    Named("intercept") = 1);
  int nspde = spde["n.spde"];
  Rcpp::List effect_list = Rcpp::List::create(
    Rcpp::List::create(Named("i") = Rcpp::seq(1,nspde)),
    Rcpp::List::create(Named("Intercept") = 1,
                       Named("coords") = dfout)
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
  Rcpp::List par_list = Rcpp::List::create(
    Named("initial") = 1.0,
    Named("fixed") = Rcpp::LogicalVector(TRUE)
  );
  Rcpp::List prec_list = Rcpp::List::create(
    Named("prec") = par_list
  );
  Rcpp::List hyper_list = Rcpp::List::create(
    Named("hyper") = prec_list
  );
  Rcpp::List inla_results = inla(Rcpp::_["formula"] = form(formula),
                                 Rcpp::_["control.family"] = hyper_list,
                                 Rcpp::_["control.compute"] = Rcpp::List::create(Named("dic") = Rcpp::LogicalVector(TRUE)),
                                 Rcpp::_["data"] = inlaStackData(stkDat, Rcpp::_["spde"] = spde),
                                 Rcpp::_["control.predictor"] = controlpred_list);
  Rcpp::NumericVector v = inla_results["mlik"]; //storing the vector of mliks
  //Rcout << v[0] << std::endl;
  Function maternmat("matern.mat");
  Rcpp::NumericMatrix mat = maternmat(Rcpp::_["nu"] = 1,
                                      Rcpp::_["kappa"] = 1,
                                      Rcpp::_["coords"] = dfout,
                                      Rcpp::_["sigma2e"] = 1,
                                      Rcpp::_["sigma2x"] = 1);
  return inla_results["mlik"];
}

