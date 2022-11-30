#ifndef MYFITCONGRAPHC_ALT_H // include guard
#define MYFITCONGRAPHC_ALT_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

void myfitConGraphC_alt(const arma::mat & amat,
                        const arma::mat & S,
                        arma::mat & Shat,
                        int & iter,
                        double tol = 0.000001);

#endif /* MYFITCONGRAPHC_ALT_H */
