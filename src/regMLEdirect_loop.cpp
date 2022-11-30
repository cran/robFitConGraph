#include <RcppArmadillo.h>
#include "myFitConGraphC_alt.h"
#include "CapplyC.h"
#include "mahalanobis_fast_aux.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List regMLEdirect_loop(arma::mat & X, 
                             arma::mat amat,
                             arma::mat S_old, 
                             arma::mat & S_new, 
                             arma::vec & mu_old, 
                             arma::vec & mu_new,
                             const double df,
                             const double tol = 0.000001){
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec ss(n);
  arma::vec us(n);
  arma::mat Shat(p,p);
  int FCG_iter = 0; 
  int tot_FCG_iter = 0; 
  double ips_it = 0;
  int em_it = 0;
  amat.diag().zeros(); // Set diagonal of amat to zero
  
  while(!arma::approx_equal(S_new, S_old, "absdiff", tol)){
    S_old = S_new;
    mu_old = mu_new;
    ss = mahalanobis_fast_aux(X, mu_old, S_old);
    // u <- function(s,p,d){return((d+p)/(d+s))}
    // us <-  sapply(X=ss,FUN=u,p=p,d=df) # u instead of uTyler
    for(int ii = 0; ii < n; ii++){
      us(ii) = (df+p)/(df+ss(ii));
    }
    mu_new = arma::vectorise(arma::trans(us)*X / arma::as_scalar(arma::sum(us)));
    S_new = CapplyC(us, X, mu_new);
    myfitConGraphC_alt(amat, S_new, Shat, FCG_iter, tol);
    tot_FCG_iter += FCG_iter;
    S_new = Shat;
    ++em_it;
  }
  ips_it = static_cast<double>(tot_FCG_iter) / static_cast<double>(em_it);
  
  return Rcpp::List::create(Rcpp::Named("S_new") = S_new,
                            Rcpp::Named("mu") = mu_new,
                            Rcpp::Named("em_it") = em_it,
                            Rcpp::Named("ips_it") = ips_it);
}
