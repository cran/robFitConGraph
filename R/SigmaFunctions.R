#------------------------------------------------
#' The proportionality constant \eqn{\eta} of the t-MLE for scatter
#'
#' @param df_data A positive real number or \code{Inf}. The degrees of freedom of the
#'   data-generationg t-distribution. \code{Inf} means normal distribution.
#'
#' @param df_est A positive real number or \code{Inf}. The degrees of freedom of the
#'   t-distribution the M-estimator is derived from. \code{Inf} is the usual sample covariance
#'   (i.e. the MLE under normality).
#'
#' @param p An integer, at least 2.
#'
#' @return A real value. Returns the constant \eqn{\eta} (cf. Vogel and Tyler 2014, p. 870, Example 2).
#'    This first appeared in Tyler (1982, p. 432, Example 3) as \eqn{\sigma^{-1}}.
#'    It is also stated in Tyler (1983, p. 418) as \eqn{\sigma^{-1}_{u,g}}.
#'
#'
#' @details Let \eqn{X_1,...,X_n} be an i.i.d.\ sample from \eqn{t_{\nu,p}(\mu, S)}, i.e.,
#' a p-variate t-distribution with \eqn{\nu} degrees of freedom, location parameter \eqn{\mu}
#' and shape matrix \eqn{S}. The limit case \eqn{\nu=\infty} is allowed, where \eqn{t_{\infty,p}(\mu,S)} is
#' \eqn{N_p(\mu,S)}.
#' \cr
#' 
#' Let \eqn{\hat{S}_n} be the \eqn{t_m} MLE for scatter. Also here, \eqn{m=\infty} is allowed:
#' This is the sample covariance matrix.
#' If \eqn{\hat{S}_n} is applied to \eqn{X_1,...,X_n}, then, as \eqn{n \to \infty},
#' \eqn{\hat{S}_n} converges in probability to \eqn{\eta S}.
#' The function \code{find_eta()} returns the proportionality constant \eqn{\eta}
#' for inputs \eqn{\nu, m} and \eqn{p}.
#' (Note: if \eqn{\nu \neq m}, \eqn{\hat{S}_n} is technically not an MLE, but an M-estimator.)
#' \cr
#' \cr
#' Some specific values:
#' \itemize{
#' \item
#'   If \eqn{\nu = m} (also for \eqn{\infty = \infty}), then \eqn{\eta = 1}
#'   (i.e., the MLE at the corresponding population
#'   distribution consistently estimates its population value).
#' \item
#'   If \eqn{m = \infty} and \eqn{2 < \nu < \infty}, then \eqn{\eta = \nu/(\nu-2)}.
#' \item
#'   If \eqn{m = \infty} and \eqn{\nu <= 2}, then \eqn{\eta = \infty}.
#'   Precisely: \eqn{\hat{S}_n} does not
#'   converge in this case.
#' }
#'
#' The general expressions: \eqn{\eta} is the solution to
#' \cr
#' \eqn{F(\eta) = E(\phi(R/\eta)) - p = 0},
#' \cr
#' where \eqn{\phi(y) = y(m+p)/(m+y)} and \eqn{R = (X - \mu)^\top S^{-1} (X-\mu)} for
#' \eqn{X \sim t_{\nu,p}(\mu,S)}.
#' For the integral, \code{stats::integrate} is used, for finding the root the
#' function \code{stats::uniroot}.
#'
#' In general, \eqn{\nu} (\code{df_data}) and \eqn{m} (\code{df_est}) can take on any
#' positive value, including \eqn{\infty}. The function works well for \eqn{p <= 100} and
#' \eqn{\nu >= 1}. For larger values of \eqn{p}, setting \eqn{\eta = 1} provides a good approximation
#' (unless df_data is very small and df_est is rather large).
#' For smaller values of \eqn{\nu}, try the following potential remedies:
#' \itemize{
#' \item Re-consider if \eqn{\nu < 1} is really necessary. This is VERY heavy-tailed.
#' \item Adaptation of the search interval.
#' \item By a suitable substitution, transform the integral to numerically more stable one (bounded support).
#' }
#'
#' \eqn{F(\eta)} is a decreasing function.
#' \itemize{
#' \item The larger df_est, the larger \eqn{\eta}.
#' \item The larger df_data, the smaller \eqn{\eta}.
#' \item The larger \eqn{p}, the closer \eqn{\eta} is to 1. }
#'
#' @author Daniel Vogel
#'
#' @references
#'   Tyler, D. E. (1982): Radial estimates and the test for sphericity,
#'   \emph{Biometrika}, 69, 2, pp. 429-36
#'   \cr\cr
#'   Tyler, D. E. (1983): Robustness and efficiency properties of scatter matrices,
#'   \emph{Biometrika}, 70, 2, pp. 411-20
#'   \cr\cr
#'   Vogel, D., Tyler, D. E. (2014): Robust estimators
#'   for nondecomposable elliptical graphical models, \emph{Biometrika}, 101, 865-882
#'
#'
#' @examples
#' find_eta(df_data = Inf, df_est = 3,   p = 10)
#' find_eta(df_data = 4.5, df_est = 4.5, p = 2)
#' @export


find_eta <- function(df_data, df_est, p){
  nu <- df_data
  m <- df_est

  if (nu <= 0) stop("df_data must be positive real or Inf")
  if (m <= 0) stop("df_est must be positive real or Inf")

  if (nu==Inf){  # normal data
    if (m==Inf)  # sample covariance at normal data
      return(1)
    else{        # t_m MLE at normal data
      integrand <- function(y,eta){ return( (m+p) / (m*eta + y) * y^(p/2) * exp(-y/2)) }
      K <- (0.5)^(p/2) / gamma(p/2)
    }
  }else{         # t_nu data
    if (m==Inf){ # the sample covaiance at t_nu data
      if (nu > 2)
        return(nu/(nu-2))
      else
        return(Inf)
    }
    else if(m==nu) {
      return(1)  # t_nu MLE at t_nu data
    }
    else{        # t_m MLE at t_nu data
      integrand <- function(y,eta){ return( (m+p) * y^(p/2) * (1 + y/nu)^(-(nu+p)/2) / (m*eta + y) ) }
      K <- c.nu.p(nu, p) * pi^(p/2) / gamma(p/2)
    }
  }

  gesInt <- function(eta){
    integr <- function(x)(integrand(x,eta))
    int <- stats::integrate(f=integr,lower=0,upper=Inf,
                            subdivisions=20000,
                            abs.tol=1e-7,
                            stop.on.error=FALSE)
                            # last option fairly important: integral may be "too large" for some
                            # values of eta; but still root (where eta near zero) can be found
    return(K*int$value-p)
  }

  eta <- stats::uniroot(f=gesInt,interval=c(0.00001,4),
                        tol = .Machine$double.eps^0.5,
                        extendInt="downX")
                        # the last option is important:
                        # eta may lie beyond the upper interval limit
                        # this fixes this nicely
  return(eta$root)
}







#------------------------------------------------
#' The asymptotic efficiency constant \eqn{\sigma_1} of the t-MLE for scatter
#'
#' @param df_data A positive real number or \code{Inf}. The degrees of freedom of the
#'   data-generationg t-distribution. \code{Inf} means normal distribution.
#'
#' @param df_est A non-negative real number or \code{Inf}. The degrees of freedom of the
#'   t-distribution the M-estimator is derived from.
#'   \code{Inf} is the usual sample covariance, \code{0} is Tyler's M-estimator.
#'
#' @param p An integer, at least 2.
#'
#' @return A real value. Returns the constant \eqn{\sigma_1} (cf. Vogel and Tyler 2014, p. 870, Example 2).
#'   This first appeared in Tyler (1982, p. 432, Example 3).
#'
#' @details 
#' Let \eqn{X_1,...,X_n} be an i.i.d. sample from \eqn{t_{\nu,p}(\mu, S)}, i.e.,
#' a p-variate t-distribution with \eqn{\nu} degrees of freedom, location parameter \eqn{\mu}
#' and shape matrix \eqn{S}. The limit case \eqn{\nu = \infty} is allowed, where \eqn{t_{\infty,p}(\mu,S)} is
#' \eqn{N_p(\mu,S)}.
#' \cr
#' \cr
#' Let \eqn{\hat{S}_n} be the \eqn{t_m} MLE for scatter.
#' Also here, \eqn{m=\infty} is allowed: This is the sample covariance matrix.
#' If \eqn{\hat{S}_n} is applied to \eqn{X_1,...,X_n}, then, as \eqn{n \to \infty},
#' \eqn{\hat{S}_n} converges in probability to \eqn{\eta S}.
#' The function \code{find_sigma1()} returns a scalar appearing in the asymptotic
#' covariance matrix of \eqn{\hat{S}_n}.
#' \cr
#' \cr
#' The scalar  \eqn{\sigma_1} is defined as 
#' \deqn{ \sigma_1  =  \frac{(p+2)^2 \gamma_1}{(2\gamma_2 + p)^2}, } 
#' where
#' \deqn{\gamma_1 = \frac{E\{\phi^2(R/\eta)\}}{p(p+2)} \quad \mbox{ and } \quad
#' \gamma_2 = \frac{1}{p} E\left\{\frac{R}{\eta}\phi'\left(\frac{R}{\eta}\right)\right\},}
#' furthermore
#' \eqn{\phi(y) = y(m+p)/(m+y)} and \eqn{R = (X - \mu)^\top S^{-1} (X-\mu)} for
#' \eqn{X \sim t_{\nu,p}(\mu,S)}, and \eqn{\eta} is defined in the help page of 
#' \code{\link{find_eta}}.
#' \cr
#' \cr
#' A noteworthy difference between \code{find_sigma1} and
#' \code{\link[robFitConGraph]{find_eta}} is that the argument \code{df_est} may be
#' \code{0} for \code{find_sigma1}, but must strictly positive for \code{find_eta}.
#' For both functions, \code{df_data} must be strictly positive. There is no such thing
#' as a t-distribution with zero degrees of freedom. There is such a thing as a
#' t-MLE with zero degrees of freedom: the Tyler estimator. Its \eqn{\sigma_1} value is
#' \eqn{1 + 2/p} regardless of the underlying elliptical distribution. However, since
#' the Tyler estimator provides shape information only, but none on scale,
#' \eqn{\eta} is irrelevant in this case.
#'
#' @author Daniel Vogel
#'
#' @references Vogel, D., Tyler, D. E. (2014): Robust estimators
#'   for nondecomposable elliptical graphical models, \emph{Biometrika}, 101, 865-882
#'   \cr\cr
#'   Tyler, D. E. (1982): Radial estimates and the test for sphericity,
#'   \emph{Biometrika}, 69, 2, pp. 429-36
#   \cr\cr
#   Tyler, D. E. (1983): Robustness and efficiency properties of scatter matrices,
#   \emph{Biometrika}, 70, 2, pp. 411-20
#'
#' @examples
#' find_sigma1(df_data = Inf, df_est = 3,   p = 10)
#' find_sigma1(df_data = 4.5, df_est = 4.5, p = 2)
#'
#' @export
find_sigma1 <- function(df_data, df_est, p){
  if (df_data <= 0)
    stop("df_data must be positive. df_est = 0 is allowed (Tyler estimator)")
  if (df_est == 0)
    sigma1 <- 1 + 2/p
  else{
    eta    <- find_eta(    df_data=df_data, df_est=df_est, p=p)
    gamma1 <- find_gamma1( df_data=df_data, df_est=df_est, p=p, eta=eta)
    gamma2 <- find_gamma2( df_data=df_data, df_est=df_est, p=p, eta=eta)
    sigma1 <- sigma1_from_gammas(gamma1 = gamma1, gamma2 = gamma2, p=p)
  }
  return(sigma1)
}







#-------------------------------------------------------------------------------------
# further auxiliary functions below


#--------------------------------
# the constant for the t-density
# nu: the degrees of freedom, real, positive
# p:  dimension, positive integer

c.nu.p <- function(nu,p){return(gamma((nu+p)/2)/ (sqrt(pi*nu))^p /gamma(nu/2))}



#--------------------------------
# the constant gamma1 (cf. Vogel and Tyler 2014, p. 870, Example 2)
# df_data, df_est: positive or Inf (0 not allowed)

find_gamma1 <- function(df_data, df_est, p, eta){
  nu <- df_data
  m <- df_est

  if (nu==Inf){   # normal data
    if (m==Inf){  # sample covariance at normal data
      integrand <- function(y,eta){ return( eta^(-2)  *  y^(p/2 + 1)  *  exp(-y/2)  ) }
      K <- 1 / ( 2^(p/2) * gamma(p/2) * p * (p+2))
    }
    else{         # t_m MLE at normal data
      integrand <- function(y,eta){ return( ((m+p)/(m*eta + y))^2  *  y^(p/2 + 1)  *   exp(-y/2) ) }
      K <-  1 / ( 2^(p/2) * gamma(p/2) * p * (p+2))
    }
  }else{          # t_nu data
    if (m==Inf){  # sample covariance at t_nu data
      if (nu < 5)
        return(NA)
      else{
        integrand <- function(y,eta){ return( eta^(-2)  *  y^(p/2 + 1)  *  (1 + y/nu)^(-(nu+p)/2)  ) }
        K <- c.nu.p(nu, p) * pi^(p/2) / gamma(p/2) / p / (p+2)
      }
    }
    else{        # t_m MLE at t_nu data
      integrand <- function(y,eta){ return( ((m+p)/(m*eta + y))^2  *  y^(p/2 + 1)  *   (1 + y/nu)^(-(nu+p)/2)  ) }
      K <- c.nu.p(nu, p) * pi^(p/2) / gamma(p/2) / p / (p+2)
    }
  }

  integr <- function(x)(integrand(x,eta))
  int <- stats::integrate(f=integr,lower=0,upper=Inf,subdivisions=20000,abs.tol=1e-7)
  return(K*int$value)
}





#--------------------------------
# the constant gamma2 (cf. Vogel and Tyler 2014, p. 870, Example 2)
# df_data, df_est: positive or Inf (0 not allowed)

find_gamma2 <- function(df_data, df_est, p, eta){
  nu <- df_data
  m <- df_est

  if (nu==Inf){    # normal data
    if (m==Inf){   # sample covariance at normal data
      integrand <- function(y,eta){ return( eta^(-1) * y^(p/2)  *  exp(-y/2)  ) }
      K <- 1 / ( 2^(p/2) * gamma(p/2) * p)
    }
    else{          # t_m MLE at normal data
      integrand <- function(y,eta){ return(  eta * m * (m+p) / (m*eta + y)^2 *  y^(p/2)  *   exp(-y/2) ) }
      K <-  1 / ( 2^(p/2) * gamma(p/2) * p )
    }
  }else{           # t_nu data
    if (m==Inf){   # sample covariance at t_nu data
      if (nu < 5)
        return(NA)
      else {
        integrand <- function(y,eta){ return(  eta^(-1)  *  y^(p/2)  *   (1 + y/nu)^(-(nu+p)/2)  ) }
        K <- c.nu.p(nu, p) * pi^(p/2) / gamma(p/2) / p
      }
    }
    else{          # t_m MLE at t_nu data
      integrand <- function(y,eta){ return(  eta * m * (m+p) / (m*eta + y)^2  *  y^(p/2)  *   (1 + y/nu)^(-(nu+p)/2)  ) }
      K <- c.nu.p(nu, p) * pi^(p/2) / gamma(p/2) / p
    }
  }

  integr <- function(x)(integrand(x,eta))
  int <- stats::integrate(f=integr,lower=0,upper=Inf,subdivisions=20000,abs.tol=1e-7)
  return(K*int$value)
}


#--------------------------------
# computes sigma1 from gamma1, gamma2 and p
sigma1_from_gammas <- function(gamma1,gamma2,p){
  if (is.na(gamma1) || is.na(gamma2))
    return(Inf)
  else
    return( (p+2)^2 * gamma1 / (2*gamma2 + p)^2 )
}



