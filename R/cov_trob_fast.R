cov_trob_fast <- function(x, wt = rep(1,n), cor = FALSE, center = TRUE, nu = 5,
                        maxit = 25, tol = 0.01)
{
  x <- as.matrix(x)
  n = nrow(x)
  p = ncol(x)
  dn <- colnames(x)
  center <- as.matrix(center)
  
  if(any(is.na(x)) || any(is.infinite(x)))
    stop("missing or infinite values in 'x'")
  
  if(missing(wt)){
    L <- aux_My_Cov_Trob(x = x,
                         center = center,
                         wtMissing = TRUE,
                         wt = wt,
                         cor = cor,
                         nu = nu,
                         maxit = maxit,
                         tol = tol)
    ans <- list(cov = L$cov,
                center = as.vector(L$loc),
                n.obs = n)
  }
  else
  {
    L <- aux_My_Cov_Trob(x = x,
                         center = center,
                         wtMissing = FALSE,
                         wt = wt,
                         cor = cor,
                         nu = nu,
                         maxit = maxit,
                         tol = tol)
    ans <- list(cov = L$cov,
                center = as.vector(L$loc),
                wt = L$wt,
                n.obs = n)
  }
  
  switch(L$err,
         stop("length of 'wt' must equal number of observations"),
         stop("negative weights not allowed"),
         stop("no positive weights present"),
         stop("'center' is not the right length"),
         stop("'center' must be a logical value or numeric vector"),
         warning("Probable convergence failure"))
  
  if(length(dn)) {
    dimnames(ans$cov) <- list(dn, dn)
    names(ans$center) <- dn
  }
  
  if(cor) {
    sd <- sqrt(diag(ans$cov))
    cor <- (ans$cov/sd)/rep(sd, rep.int(p, p))
    if(length(dn)) dimnames(cor) <- list(dn, dn)
    ans <- c(ans, list(cor = cor))
  }
  ans$call <- match.call()
  ans$iter <- L$endit
  ans
}
