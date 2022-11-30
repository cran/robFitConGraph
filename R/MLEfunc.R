# April 2022
# Daniel Vogel

#' Graph-constrained robust scatter estimation.
#'
#' The function computes a robust estimate of a scatter matrix subject to
#' zero-constraints in its inverse. The methodology is described in Vogel &
#' Tyler (2014).
#'
#' @param X A data matrix with \eqn{n} rows and \eqn{p} columns, representing
#'   \eqn{n} observations and \eqn{p} variables. Elements of \code{X} must be
#'   numeric and \eqn{n} must be at least \eqn{p+1}. Rows containing NA are omitted.
#'
#' @param amat A \eqn{p} times \eqn{p} matrix representing the adjacency matrix
#'   of a graphical model. \code{amat} must be symmetric with numerical entries
#'   \code{0} or \code{1}. Alternatively, a Boolean matrix. The entries on the diagonal 
#'   are irrelevant. Must not contain any NAs.
#'
#' @param df A positive real number or \code{Inf}. The degrees of freedom of the t-distribution
#'   used (see Details below). The value \code{df = Inf} corresponds to the sample 
#'   covariance matrix. The value \code{df = 0} is not allowed as Tyler's M-estimator is 
#'   currently not implemented. Default value is \code{3}.
#'   
#' @param tol tolerance for numerical convergence. Iteration stops if the
#'   maximal element-wise difference between two successive matrices is less
#'   than \code{tol}. Must be at least \code{10e-14}. Default is \code{1e-5}.
#'
#' @param sigma1 numerical. A positive real value. This value is needed for the p-value of the model fit.
#'   If missing, it is computed by \code{find_sigma1(df_data = Inf, df_est = df, p = ncol(X))}.
#'   See Details below and the documentation of the function \code{\link[robFitConGraph]{find_sigma1}}.
#'
#' @param plug_in logical. The function offers two types of estimates: the
#'   plug-in M-estimator and the direct M-estimator. If \code{plug_in} is
#'   \code{TRUE}, the plug-in estimate is computed. If \code{FALSE}, the direct
#'   M-estimator is computed. The plug-in estimator is faster, but has higher variance.
#'   Default is \code{TRUE}. Ignored if \code{df == Inf}.
#'
#' @param direct logical. If \code{TRUE}, the direct estimate is computed,
#'   otherwise the plug-in estimate. Default is \code{FALSE}. In case of
#'   conflicting specifications of \code{plug-in} and \code{direct}, \code{plug_in} overrides
#'   \code{direct}. Ignored if \code{df == Inf}.
#'
#' 
#' @return List with 8 elements:
#'   \item{\code{Shat}}{\eqn{p \times p}{p times p} symmetric scatter matrix estimate.}
#'   \item{\code{mu}}{numerical \eqn{p}-vector, the location estimate
#'     In case of \code{df == Inf}, this is the sample mean.}
#'   \item{\code{pval}}{numerical. The p-value of the model fitted. \eqn{D_n/\sigma_1}
#'   is compared to the \eqn{\chi_r^2} distribution, see Details.}
#'   \item{\code{dev}}{numerical. The deviance test statistic \eqn{D_n}, see Details.}
#'   \item{\code{missing_edges}}{integer. The number of missing edges \eqn{r} as indicated
#'   by the argument \code{amat}.}
#'   \item{\code{sigma1}}{numerical. Either the argument \code{sigma1} (if provided), or the
#'   return value of
#'   \code{find_sigma1(df_data = Inf, df_est = df, p = ncol(X))}}
#'   \item{\code{em_it}}{integer. Number of iterations of the t-MLE computation.}
#'   \item{\code{ips_it}}{integer. In the case of the plug-in estimate, this is
#'        the number of iterations of the Gaussian graphical model fitting procedure
#'   (Algorithm 17.1) in Hastie et al 2004). In the case of the direct estimate,
#'   the Gaussian graphical model fitting is executed \code{em_it} times, and \code{ips_it} returns
#'   the average number of iterations.}
#'
#'    The last five return values are 
#'    mainly for traceability purposes.
#'    Particularly \code{dev}, \code{missing_edges}, and \code{sigma1} 
#'    are easily obtained by the inputs or outputs of \code{robFitConGraph}.
#'    They are the ingredients needed to compute the p-value:
#'
#'    \code{
#'      pval <- pchisq(q = dev/sigma1, df = missing_edges, lower.tail = FALSE)
#'    }
#'
#'
#'
#' @details Two types of graph-constrained M-estimates of scatter based 
#'   on the elliptical t-distribution are implemented: the direct estimate and the
#'   plug-in estimate. The direct estimate is referred to as graphical
#'   M-estimator in Vogel & Tyler (2014).
#'
#'   The plug-in estimate is two algorithms performed sequentially: First an
#'   unconstrained t-maximum likelihood estimate of scatter is computed (the
#'   same as \code{\link[MASS]{cov.trob}} from \code{MASS}). This is then
#'   plugged into the Gaussian graphical model fitting routine (the same as
#'   \code{\link[ggm]{fitConGraph}} from \code{ggm}). Specifically
#'    Algorithm 17.1 from Hastie, Tibshirani, Friedman (2009) is used.
#'
#'   The direct estimate is the actual maximum-likelihood estimator within the
#'   elliptical graphical model based on the elliptical t-distribution. The
#'   algorithm is an iteratively-reweighted least-squares algorithm, where the
#'   Gaussian graphical model fitting procedure is nested into the t-estimation
#'   iteration. The direct estimate therefore takes longer to compute, but the
#'   estimator has a better statistical efficiency for small sample sizes. Both
#'   estimators are asymptotically equivalent. The estimates tend to be very
#'   close to each other for large sample sizes.
#'
#'   \bold{The deviance test statistic}
#'
#'   The robustified deviance test statistic (short: deviance) for testing 
#'   a graphical model \eqn{G} is
#'   \deqn{
#'       D_n = n ( \log \det(\hat{S}_n(G)) - \log \det(\hat{S}_n) ),  }{%
#'       D_n = n ( log det(S_n(G)) - log det(S_n) ),  }
#'   where \eqn{\hat{S}_n(G)}{S_n(G)} is a graph-constrained scatter estimator and \eqn{\hat{S}_n}{S_n} the
#'   corresponding unconstrained scatter estimator, see Vogel & Tyler (2014, p. 866 bottom).
#'   Under the graphical model \eqn{G}, the deviance \eqn{D_n} converges to
#'   \eqn{\sigma_1 \chi_r^2}, where \eqn{r} is the number of missing edges in \eqn{G}. The
#'   constant \eqn{\sigma_1} depends on the scatter estimate \eqn{\hat{S}_n}{S_n}, the dimension \eqn{p},
#'   and the population distribution. It can either be provided directly (as \code{sigma1})
#'   or, if missing, is computed by \code{\link{find_sigma1}}, assuming a Gaussian population
#'   distribution.
#'
#'   Although \code{robFitConGraph} combines the functionality of
#'   \code{\link[ggm]{fitConGraph}} and \code{\link[MASS]{cov.trob}} and
#'   contains both as special cases, it uses neither. The
#'   algorithms are implemented in C++. 
#'
#'   Input and output of \code{robFitConGraph} are intended to resemble
#'   \code{\link[ggm]{fitConGraph}} from the package \code{ggm}.
#'   Some notable differences:
#'
#'   \itemize{
#'   \item \code{\link[ggm]{fitConGraph}} takes as input the
#'      unconstrained covariance matrix; \code{robFitConGraph} takes the actual data.
#'
#'   \item \code{\link[ggm]{fitConGraph}} alternatively offers to  specify the
#'      graphical model as list of cliques; \code{robFitConGraph} only takes
#'      the adjacency matrix \code{amat}.
#'
#'    \item \code{\link[ggm]{fitConGraph}} returns a value called the \emph{degrees of freedom} as \code{df}. 
#'      These degrees of freedom mean the number of missing edges, which appear 
#'      as the degrees of freedom of the chi-square distribution of the deviance test.
#'      However, these degrees of freedom are completely unrelated to the
#'      the argument \code{df_est} of \code{robFitConGraph}, which refers to the
#'      degrees of freedom of the t-distribution defining the scatter estimate.
#'      \code{robFitConGraph} returns the number of missing edges as \code{missing_edges}. 
#'    }
#'
#' @author Stuart Watt, Daniel Vogel
#'
#' @references Vogel, D., Tyler, D. E. (2014): Robust estimators
#'   for nondecomposable elliptical graphical models, \emph{Biometrika}, 101, 865-882\cr\cr
#'   Hastie, T., Tibshirani, R. and Friedman, J. (2004). \emph{The elements of
#'   statistical learning}. New York: Springer.
#'
#' @seealso \code{\link[ggm]{fitConGraph}} from package
#'   \code{ggm} for non-robust graph-constrained covariance estimation
#'   \cr
#'   \code{\link[MASS]{cov.trob}} from package \code{MASS} for unconstrained
#'   \code{p} times \code{p} t-MLE scatter matrix
#'

#' @examples
#' # --- Example 1: anxieties data ---
#'
#' # The data set 'anxieties' contains 8 anxiety variables of 82 observations. We test the
#' # null hypothesis that, given Generalized Anxiety (GAD), Music performance anxiety (MPA)
#' # is conditionally independent of all remaining 6 variables:
#'
#' # - build the corresponding graphical model -
#' data(anxieties)
#' amat = matrix(1,ncol=ncol(anxieties),nrow=ncol(anxieties))
#' colnames(amat) <- colnames(anxieties)
#' rownames(amat) <- colnames(amat)
#' amat["MPA",c("SAD","PD","AG","SP","SEP","ILL")] <- 0
#' amat[c("SAD","PD","AG","SP","SEP","ILL"),"MPA"] <- 0
#' amat
#' #     MPA GAD SAD PD AG SP SEP ILL
#' # MPA   1   1   0  0  0  0   0   0
#' # GAD   1   1   1  1  1  1   1   1
#' # SAD   0   1   1  1  1  1   1   1
#' # PD    0   1   1  1  1  1   1   1
#' # AG    0   1   1  1  1  1   1   1
#' # SP    0   1   1  1  1  1   1   1
#' # SEP   0   1   1  1  1  1   1   1
#' # ILL   0   1   1  1  1  1   1   1
#'
#' robFitConGraph(X = anxieties, amat = amat)$pval
#' # [1] 0.881159
#' # This provides no evidence against the null hypothesis.
#'
#' # the non-robust, sample-covariance-based model fit:
#' robFitConGraph(X = anxieties, amat = amat, df = Inf)$pval
#' # [1] 0.6310895
#' # ...provides no evidence against the null hypothesis either.
#' 
#' # The latter is also obtained as:
#' # pchisq(q = ggm::fitConGraph(S = cov(anxieties), amat = amat, 
#' #                             n = nrow(anxieties))$dev, 
#' #        df = 6, lower.tail = FALSE)
#'
#'
#'
#' # --- Example 2: simulated data - the chordless 7-cycle ---
#'
#' # - build the corresponding graphical model -
#' chordless_p_cycle <- function(rho, p){
#'   M <- diag(1,p)
#'   for (i in 1:(p-1)) M[i,i+1] <- M[i+1,i] <- -rho
#'   M[1,p] <- M[p,1] <- -rho
#'   return(M)
#' }
#' p <- 7                             # number of variables
#' rho <- 0.4                         # partial correlation
#' PCM <- chordless_p_cycle(rho, p)    # partial correlation matrix
#' SM <- cov2cor(solve(PCM))          # shape matrix (i.e covariance matrix up to scale)
#' model <- abs(sign(PCM))            # adjacency matrix of the chordless-7-cycle
#' # > model
#' #      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#' # [1,]    1    1    0    0    0    0    1
#' # [2,]    1    1    1    0    0    0    0
#' # [3,]    0    1    1    1    0    0    0
#' # [4,]    0    0    1    1    1    0    0
#' # [5,]    0    0    0    1    1    1    0
#' # [6,]    0    0    0    0    1    1    1
#' # [7,]    1    0    0    0    0    1    1
#'
#' # This is the cordless-7-cycle (p.872 Figure 1 (a) in Vogel & Tyler, 2014).
#' # All non-zero partial correlations are 0.4.
#' # The true covariance is (up to scale) 'SM'. This matrix is constructed such
#' # that it has zero entries in its inverse as specified by 'model'.
#'
#'
#' # - generate data from the graphical model (elliptical t3 distribution) -
#' n <- 50            # number of observations
#' df_data <- 3       # degrees of freedom
#' library(mvtnorm)   # for rmvt function
#' set.seed(918273)   # for reproducability
#' X <- rmvt(n = n, sigma = SM, df = df_data)
#'
#' # X contains a data set of size n = 50 and dimension p = 7, sampled from the
#' # elliptical t-distribution with 3 degrees of freedom and shape matrix 'SM'
#'
#'
#' # - comparing estimates -
#'
#' # the true correlation matrix:
#' round(cov2cor(SM), d = 2)
#'
#' # Pearson correlations:
#' round(cov2cor(cov(X)), d = 2)
#'
#' # robust correlations based on the direct graph-constrained t-MLE:
#' S1 <- robFitConGraph(X, amat = model, df = df_data, direct = TRUE)$Shat
#' round(cov2cor(S1), d = 2)
#'
#' # robust correlations based on the plug_in graph-constrained t-MLE:
#' S2 <- robFitConGraph(X, amat = model, df = df_data, plug_in = TRUE)$Shat
#' round(cov2cor(S2), d = 2)
#'
#' # The correlation estimates based on S1 and S2 are close to the true
#' # correlations, whereas the sample correlations differ strongly.
#' # Note: sample correlations are not asymptotically normal at the t3 distribution.
#' @export

robFitConGraph <- function(X, amat, df = 3, tol = 1e-5, sigma1,
                              plug_in = TRUE, direct = FALSE){

  # =========================================================================================================
  #
  #   Initial checks & preliminaries 
  #   - stop checks (errors) and messages
  #   - na.omit on X
  #   - p and n defined
  #   - amat diagonal set to 0
  #
  # =========================================================================================================

  #- X check -
  if(missing(X)) stop('data matrix X required')
  if(is.numeric(X) == FALSE) stop('data matrix X must be numeric data type')
  if (any(is.na(X))){
    X <- stats::na.omit(X)
    message('There are NAs in X. Incomplete lines are removed.')
  }

  #- amat check -
  if(missing(amat)) stop('amat matrix needs to be specified')
  if(isSymmetric(amat) == FALSE) stop('amat is not symmetric')
  if(length(amat[col(amat)!=row(amat)][!(amat[col(amat)!=row(amat)] %in% c(0,1))]) != 0) stop('off-diagonal amat entries must be 0 or 1')
  if (!mode(amat) %in% c("numeric","integer","logical","double")) stop('amat must be of mode numeric or logical.')
  if (any(is.na(amat))) stop('amat contains NAs.')
  
  #- fixing amat to a numeric 0-1-matrix with 0 on the diagonal
  # important for computing the number of missing edges correctly below
  amat <- amat*1.0
  diag(amat) <- 0


  #- dimension check -
  if((ncol(X) != ncol(amat))) stop('number of columns of X and amat must match')
  if(ncol(amat) < 2) stop('number of variables must be greater than or equal to 2')
  if(nrow(X) <= ncol(X)) stop('number of observations must be greater than number of variables')
  n <- nrow(X)
  p <- ncol(X)

  #- df check -
  if(is.numeric(df) == FALSE) stop('data input df must be double data type')
  if(df < 0){
      message('df provided is negative, but must be positive. Replaced by default df = 3.')
      df <- 3
  } else if (df == 0) {
      message("df provided is zero, but must be positive. Replaced by default df = 3.")
      df <- 3
  }

  #- tol check -
  if(is.numeric(tol) == FALSE) stop('data input tol must be double data type')
  if((tol < 10e-14)) stop('data input tol must be >= 10e-14')

  #- algorithm specification check -
  # plug_in_int is the internal plug_in switch
  # it used below to branch
  
  if (df != Inf){
  # if df == Inf, inputs plug_in and direct are completely ignored
      if (missing(plug_in)){
        if (missing(direct)) 
          plug_in_int <- TRUE
        else if (isTRUE(direct)) 
          plug_in_int <- FALSE
        else if (isFALSE(direct))
          plug_in_int <- TRUE
        else {
          message('Unrecognizable specification for direct. plug_in = TRUE is chosen.')
          plug_in_int <- TRUE
        }
      } else if (isTRUE(plug_in)){
        if (missing(direct)) 
          plug_in_int <- TRUE
        else if (isTRUE(direct)){
          message('Conflicting specifications for plug_in and direct. plug_in = TRUE is chosen.')
          plug_in_int <- TRUE
        } else if (isFALSE(direct))
          plug_in_int <- TRUE
        else {
          message('Unrecognizable specification for direct. plug_in = TRUE is chosen.')
          plug_in_int <- TRUE
        }    
      } else if (isFALSE(plug_in)){
        if (missing(direct)) 
          plug_in_int <- FALSE
        else if (isTRUE(direct))
          plug_in_int <- FALSE
        else if (isFALSE(direct)){
          message('Conflicting specifications for plug_in and direct. plug_in = FALSE is chosen.')
          plug_in_int <- FALSE
        } else {
          message('Unrecognizable specification for direct. plug_in = FALSE is chosen.')
          plug_in_int <- FALSE
        }      
      } else {  # branch: plug_in is something else than TRUE or FALSE
        if (missing(direct)){
          message('Unrecognizable specification for plug_in. plug_in = TRUE is chosen.')
          plug_in_int <- TRUE
        } else if (isTRUE(direct)){
          message('Unrecognizable specification for plug_in. plug_in = FALSE is chosen.')
          plug_in_int <- FALSE
        } else if (isFALSE(direct)){
          message('Unrecognizable specification for plug_in. plug_in = TRUE is chosen.')
          plug_in_int <- FALSE
        } else {
          message('Unrecognizable specifications for plug_in and direct. plug_in = TRUE is chosen.')
          plug_in_int <- TRUE
        }      
      }
  }
  
  #- sigma1 check -
  if (missing(sigma1)){
    sigma1 <- find_sigma1(df_est=df, df_data=Inf, p=ncol(X))
  }

  # =========================================================================================================
  #
  #   checks complete
  #
  # =========================================================================================================
  
  
  if(df == Inf){
      # =========================================================================================================
      #
      #   df = Inf
      #   - this is the Gaussian MLE (the usual non-robust) -
      #   - plug_in and direct coincide -
      #
      # =========================================================================================================

      S0 = stats::cov(X)
      L <- myfitConGraphC(S = S0, amat = amat, tol = tol)
      S_new  <- L$Shat
      mu     <- colMeans(X)
      em_it  <- 0
      ips_it <- L$iter

   } else {
     # =========================================================================================================
     #
     #   Regular t-MLE
     #    - this is t-MLE -
     #    - with degrees of freedom 0 < df < Inf -
     #    - distinction btw plug-in and direct case -
     #
     # =========================================================================================================
     
    
      if (plug_in_int){
   
          Cov <- cov_trob_fast(X, nu = df, tol = tol, maxit = 1000) 
          S0 <- Cov$cov

          L <- myfitConGraphC(S = S0, amat = amat, tol = tol)
          S_new <- L$Shat
          mu <- Cov$center
          em_it <- Cov$iter
          ips_it <- L$iter


      } else {
     
          Cov <- cov_trob_fast(X, nu = df, tol = tol, maxit = 1000) 
          S0 <- Cov$cov
          
          Laux <- myfitConGraphC(S = S0, amat = amat, tol = tol)
          S_new <- Laux$Shat
          L <- regMLEdirect_loop(X = X, amat = amat, 
                                 S_old = S0, S_new = S_new, 
                                 mu_old = Cov$center, mu_new = Cov$center, 
                                 df = df, tol = tol)
          # note: S_old is passed by value (no '&' in C++ code)
          S_new <- L$S_new
          mu <- drop(L$mu)
          em_it <- L$em_it
          ips_it <- L$ips_it

      } 
  }

  #--- final computation and return ---
  # each of the if-clauses above must have the following variables defined at the end:
  # S0, S_new, mu, em_it, ips_it

  dev <- n*(log(det(S_new)) - log(det(S0)))
  missing_edges <- (p*(p-1) - sum(amat))/2 # (p*(p-1)/2 edges at most; here assumed diagonal is zero)
  pval <- stats::pchisq(q = dev/sigma1, df = missing_edges, lower.tail = FALSE)

  return(list(Shat = S_new,
              mu   = mu,
              pval = pval,
              dev  = dev,
              missing_edges = missing_edges,
              sigma1 = sigma1,
              em_it  = em_it,
              ips_it = ips_it))
}
