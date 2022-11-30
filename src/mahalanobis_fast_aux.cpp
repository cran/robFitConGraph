#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec mahalanobis_fast_aux(const arma::mat & x, const arma::vec & center, 
                               const arma::mat & cov) {
  arma::mat centered = x.each_row() - center.t();
  arma::mat aux = (centered * arma::inv(cov));
  return static_cast<arma::vec>(arma::sum(aux % centered, 1));
}
