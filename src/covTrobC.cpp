#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

Rcpp::List aux_My_Cov_Trob(arma::mat x,
                           SEXP center,
                           bool wtMissing,
                           arma::vec wt,
                           bool cor,
                           double nu,
                           int maxit,
                           double tol)
{
  int n = x.n_rows;
  int p = x.n_cols;
  
  int err = 0;
  
  arma::vec wt0(n);
  
  double wtSum = 0;
    
  if(!wtMissing)
  {
    if(static_cast<int>(wt.n_elem) != n)
    {
      return Rcpp::List::create(Rcpp::Named("err") = 1);
    }
    
    int nPositive = 0;
    
    for(int i = 0; i < n; ++i)
    {
      if(wt.at(i) < 0)
      {
        return Rcpp::List::create(Rcpp::Named("err") = 2);
      }
      if(wt.at(i) > 0)
      {
        ++nPositive;
      }
    }
    if(nPositive == 0)
    {
      return Rcpp::List::create(Rcpp::Named("err") = 3);
    }
    
    // copy wt to wt0;
    for(int i = 0; i < n; ++i)
    {
      wt0.at(i) = wt.at(i);
    }
    
    arma::mat xNew(nPositive,p);
    arma::vec wtNew(nPositive);
    int dummyCount = 0;
    for(int i = 0; i < n; ++i)
    {
      if(wt.at(i) > 0)
      {
        for(int j = 0; j < p; ++j)
        {
          xNew(dummyCount,j) = x.at(i,j);
        }
        wtNew[dummyCount] = wt.at(i);
        wtSum += wt.at(i);
        ++dummyCount;
        if(dummyCount == nPositive)
        {
          break;
        }
      }
    }
    x = xNew;
    wt = wtNew;
    
    n = nPositive;
  }
  else
  {
    wtSum = n;
    wt = arma::ones(n);
  }
  
  arma::vec loc(p);
  
  for(int i = 0; i < p; ++i)
  {
    double total = 0;
    for(int j = 0; j < n; ++j)
    {
      total += wt.at(j) * x.at(j,i);
    }
    loc.at(i) = total / wtSum;
  }
  
  bool useLoc = false;
  switch(TYPEOF(center))
  {
  case REALSXP:
  {
    arma::vec tmp = (Rcpp::as<arma::vec>(center));
    if(static_cast<int>(tmp.n_elem) != p)
    {
      return Rcpp::List::create(Rcpp::Named("err") = 4);
    }
    for(int i = 0; i < p; ++i)
    {
      loc.at(i) = tmp.at(i);
    }
    break;
  }
  case LGLSXP:
  {
    bool tmp = (Rcpp::as<bool>(center));
    if(!tmp)
    {
      for(int i = 0; i < p; ++i)
      {
        loc.at(i) = 0;
      }
    }
    else
    {
      useLoc = true;
    }
    break;
  }
  default:
  {
    return Rcpp::List::create(Rcpp::Named("err") = 5);
  }
  }
  
  arma::vec w(n);
  double wTot = 0;
  for(int i = 0; i < n; ++i)
  {
    w.at(i) = wt.at(i) * (1 + p/nu);
    wTot += w.at(i);
  }
  
  int endit = 0;
  
  arma::vec w0(n);
  arma::mat X(n,p);
  
  arma::mat sXu, sXv;
  arma::vec sXd;
  
  arma::mat svd_in;
  
  double tot;
  arma::mat wX;
  arma::vec Q;
  
  for(int i = 0; i < maxit; ++i)
  {
    w0 = w;

    X = x.each_row() - loc.t();
    
    //sX
    svd_in = arma::sqrt(w / arma::sum(w)) % X.each_col();
    arma::svd_econ(sXu, sXd, sXv, svd_in);
    
    //wX & Q
    wX = X * sXv * arma::diagmat(1/sXd);
    Q = arma::vectorise(arma::square(wX) * arma::ones<arma::vec>(p));

    w = (wt * (nu + p))/(nu + Q);
    
    wTot = arma::sum(w);
    
    if(useLoc)
    {
      for(int ii = 0; ii < p; ++ii)
      {
        tot = 0;
        for(int jj = 0; jj < n; ++jj)
        {
          tot += w.at(jj) * x.at(jj, ii);
        }
        loc.at(ii) = tot / wTot;
      }
    }
    R_CheckUserInterrupt();

    // check for convergence
    bool converge = true;
    for(int ii = 0; ii < n; ++ii)
    {
      if(std::abs(w.at(ii) - w0.at(ii)) > tol)
      {
        converge = false;
        break;
      }
    }
    if(converge)
    {
      break;
    }
    
    endit = i;
  }
  
  if((endit == maxit) ||
     (std::abs(arma::mean(w)-arma::mean(wt)) > tol) ||
     ((std::abs(arma::mean(w % Q)/p - 1)) > tol))
  {
    err = 6;
  }
  
  arma::mat auxCov = arma::sqrt(w) % X.each_col();
  arma::mat cov(p,p);
  
  for(int i = 0; i < p; ++i)
  {
    for(int j = 0; j < p; ++j)
    {
      tot = 0;
      for(int k = 0; k < n; ++k)
      {
        tot += auxCov.at(k,i) * auxCov.at(k,j);
      }
      cov.at(i,j) = tot/wtSum;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("cov") = cov,
                            Rcpp::Named("loc") = loc,
                            Rcpp::Named("n,obs") = n,
                            Rcpp::Named("wt") = wt0,
                            Rcpp::Named("endit") = endit + 1,
                            Rcpp::Named("err") = err);
}
