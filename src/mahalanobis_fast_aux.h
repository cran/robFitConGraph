#ifndef MAHALANOBIS_FAST_AUX_H // include guard
#define MAHALANOBIS_FAST_AUX_H

#include <RcppArmadillo.h>
#include "mahalanobis_fast_aux.h"

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec mahalanobis_fast_aux(const arma::mat & x, const arma::vec & center, 
                               const arma::mat & cov);

#endif /* MAHALANOBIS_FAST_AUX_H */
