#ifndef CAPPLYC_H // include guard
#define CAPPLYC_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat CapplyC(const arma::vec & us, 
                  const arma::mat & X, 
                  const arma::vec & mu);
  
#endif /* CAPPLYC_H */
