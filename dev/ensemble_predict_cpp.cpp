// ============================================================================
// ensemble_predict_cpp.cpp
//
// Fast C++ kernel for BART ensemble prediction. Drop-in replacement for the
// pure-R `ensemble_predict()` in dev/rbart.R, with identical numerical output
// (verified bit-for-bit by dev/test_ensemble_predict_cpp.R).
//
// Tree representation matches dev/rbart.R exactly: each tree is an R list with
// integer/numeric vectors keyed by node id (1-based as stored, converted to
// 0-based here). Fields used at prediction:
//   is_leaf    (LogicalVector / IntegerVector — both accepted)
//   left, right     (1-based child node ids; only valid when !is_leaf)
//   split_var       (1-based covariate column index; only valid when !is_leaf)
//   split_val       (split threshold; only valid when !is_leaf)
//   mu              (leaf value; only valid when is_leaf)
//
// Build via Rcpp::sourceCpp("dev/ensemble_predict_cpp.cpp"). No external
// dependencies beyond Rcpp itself.
// ============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// Predict one tree's contribution at all rows of X.
// Pure traversal — no allocation in the inner loop besides the output vector.
// [[Rcpp::export]]
NumericVector tree_predict_cpp(List tree, NumericMatrix X) {
  // Field extraction. Accept is_leaf as either LogicalVector or IntegerVector
  // so we don't crash on R semantics quirks.
  IntegerVector left      = tree["left"];
  IntegerVector right     = tree["right"];
  IntegerVector split_var = tree["split_var"];
  NumericVector split_val = tree["split_val"];
  NumericVector mu        = tree["mu"];

  // is_leaf can be logical or integer in R. Normalize to a bool buffer once
  // to avoid per-node SEXP dispatch in the hot loop.
  SEXP is_leaf_sexp = tree["is_leaf"];
  int n_nodes = mu.size();
  std::vector<unsigned char> is_leaf(n_nodes);
  if (TYPEOF(is_leaf_sexp) == LGLSXP) {
    LogicalVector v(is_leaf_sexp);
    for (int k = 0; k < n_nodes; k++) is_leaf[k] = (v[k] == TRUE);
  } else if (TYPEOF(is_leaf_sexp) == INTSXP) {
    IntegerVector v(is_leaf_sexp);
    for (int k = 0; k < n_nodes; k++) is_leaf[k] = (v[k] != 0);
  } else {
    stop("tree$is_leaf must be logical or integer");
  }

  const int n_rows = X.nrow();
  NumericVector out(n_rows);

  // Hot loop: per row, walk the tree until reaching a leaf.
  // Stored ids are 1-based (R convention) — subtract 1 for C++ indexing.
  for (int i = 0; i < n_rows; i++) {
    int cur = 0;  // root = node id 1 (R) = index 0 (C)
    while (!is_leaf[cur]) {
      int var_idx = split_var[cur] - 1;
      if (X(i, var_idx) < split_val[cur]) {
        cur = left[cur] - 1;
      } else {
        cur = right[cur] - 1;
      }
    }
    out[i] = mu[cur];
  }
  return out;
}

// Sum tree predictions across an ensemble. Equivalent to:
//   Reduce(`+`, lapply(ensemble, function(t) tree_predict_cpp(t, X)))
// but without allocating an intermediate list of per-tree vectors.
// [[Rcpp::export]]
NumericVector ensemble_predict_cpp(List ensemble, NumericMatrix X) {
  const int n_trees = ensemble.size();
  const int n_rows = X.nrow();
  NumericVector out(n_rows);  // zero-initialized

  if (n_trees == 0) return out;

  // Reuse a buffer across trees to avoid n_trees allocations.
  for (int t = 0; t < n_trees; t++) {
    List tree = ensemble[t];

    IntegerVector left      = tree["left"];
    IntegerVector right     = tree["right"];
    IntegerVector split_var = tree["split_var"];
    NumericVector split_val = tree["split_val"];
    NumericVector mu        = tree["mu"];
    SEXP is_leaf_sexp = tree["is_leaf"];

    int n_nodes = mu.size();
    std::vector<unsigned char> is_leaf(n_nodes);
    if (TYPEOF(is_leaf_sexp) == LGLSXP) {
      LogicalVector v(is_leaf_sexp);
      for (int k = 0; k < n_nodes; k++) is_leaf[k] = (v[k] == TRUE);
    } else if (TYPEOF(is_leaf_sexp) == INTSXP) {
      IntegerVector v(is_leaf_sexp);
      for (int k = 0; k < n_nodes; k++) is_leaf[k] = (v[k] != 0);
    } else {
      stop("tree$is_leaf must be logical or integer");
    }

    for (int i = 0; i < n_rows; i++) {
      int cur = 0;
      while (!is_leaf[cur]) {
        int var_idx = split_var[cur] - 1;
        if (X(i, var_idx) < split_val[cur]) {
          cur = left[cur] - 1;
        } else {
          cur = right[cur] - 1;
        }
      }
      out[i] += mu[cur];
    }
  }
  return out;
}
