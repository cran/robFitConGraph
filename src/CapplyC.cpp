#include <RcppArmadillo.h>
#include "CapplyC.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Like CapplyB, but returns a matrix instead of a vector and all entries of output
// are divided by n.
arma::mat CapplyC(const arma::vec & us, 
                  const arma::mat & X, 
                  const arma::vec & mu){
  int n = us.n_elem;
  int p = mu.n_elem;
  int ia, ib;
  arma::mat CapplyCout(p,p,arma::fill::zeros);
  
  for(ia = 0; ia < p*p; ia++){
    for(ib = 0; ib < n; ib++){
      CapplyCout.at(ia%p,ia/p) += us.at(ib) * (X.at(ib,ia/p) - mu.at(ia/p)) * (X.at(ib,ia%p) - mu.at(ia%p));
    }
  }
  return CapplyCout / n;
}
