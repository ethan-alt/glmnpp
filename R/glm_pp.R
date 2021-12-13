


#' Power prior with fixed a0 for generalized linear models using Stan
#'
#' Sample from the posterior of a power prior (PP) generalized linear model
#' (GLM) with a fixed a0.
#'
#' @param formula         an object of class \code{\link[stats]{formula}}.
#' @param family          an object of type \code{\link[stats]{family}} giving distribution and link function
#' @param data            a \code{\link[base]{data.frame}} of current data giving all variables in \code{formula}
#' @param histdata        a \code{\link[base]{data.frame}} of historical data giving all variables in \code{formula}
#' @param a0              fixed power prior parameter
#' @param beta0           mean for initial prior on regression coefficients. Defaults to vector of 0s
#' @param Sigma0          covariance matrix for initial prior on regression coefficients. Defaults to \code{diag(100, ncol(X))}
#' @param offset          offset for current data in GLM. If \code{NULL}, no offset is utilized
#' @param offset0         offset for historical data in GLM. If \code{NULL}, no offset is utilized
#' @param disp.shape      shape parameter for inverse-gamma prior on dispersion (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param disp.scale      scale parameter for inverse-gamma prior on dispersion (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param ...             optional parameters to pass onto `rstan::sampling`
#'
#' @return an object of class [rstan::stanfit] returned by `rstan::sampling`
#' @examples
#' ## Generate current and historical data
#' set.seed(123)
#' n  = 30
#' n0 = 20
#' x  = rnorm(n)
#' x0 = rnorm(n0)
#' y  = rbinom(n = n,  size = 1, prob = binomial()$linkinv(1 - 0.5 * x) )
#' y0 = rbinom(n = n0, size = 1, prob = binomial()$linkinv(1 - 0.5 * x0) )
#' data = data.frame('y' = y, 'x' = x)
#' histdata = data.frame('y' = y0, 'x' = x0)
#'
#' ## Sample from posterior distribution of non-normalized power prior
#' fit.pp = glm_pp(y ~ x, binomial(), data, histdata, a0 = 0.5)
#' @export
glm_pp = function(
  formula, family, data, histdata, a0 = 0.5, beta0 = NULL, Sigma0 = NULL,
  offset = NULL, offset0 = NULL, disp.shape = 2.1, disp.scale = 1.1, ...
) {
  ## make sure formula is two-sided
  if(!(formula.tools::is.two.sided(formula))) {
    stop("formula must be two-sided")
  }
  ## check power prior input
  if ( a0 < 0 | a0 > 1)
    stop("Power prior parameter a0 should be between 0 and 1")
  if(any(is.na(data)))     { stop('There are missing values in data') }
  if(any(is.na(histdata))) { stop('There are missing values in histdata') }

  ## get design matrices
  X  = model.matrix(formula, data)
  X0 = model.matrix(formula, histdata)

  ## get responses
  y  = data[, all.vars(formula)[1]]
  y0 = histdata[, all.vars(formula)[1]]

  ## get mu link function as integer
  links = c('identity', 'log', 'logit', 'inverse', 'probit', 'cauchit', 'cloglog', 'sqrt', '1/mu^2')
  mu_link = which(links == family$link)[1]
  if ( length(mu_link) == 0 ) { stop(paste('Link must be one of', paste(links, collapse = ', '))) }

  ## check if offset is included
  incl_offset = 1
  if ( is.null(offset) ) {
    incl_offset = 0
    offset  = rep(0, length(y))
    offset0 = rep(0, length(y0))
  }

  ## check dimensions of offset
  if ( length(offset) != length(y) )   { stop("length of offset must match data") }
  if ( length(offset0) != length(y0) ) { stop("length of offset0 must match histdata") }

  ## default mean hyperparmeter for initial prior on beta is 0 and check dimensions
  if(is.null(beta0)) { beta0 = rep(0, ncol(X)) }
  if ( length(beta0) != ncol(X) ) { stop("specified beta0 does not match dimension of design matrix") }
  ## default covariance for initial prior on beta is diag(100)
  if(is.null(Sigma0)) { Sigma0 = diag(100, ncol(X)) }
  if( (ncol(Sigma0) != length(beta0)) | ( nrow(Sigma0) != length(beta0)) ) {
    stop('specified Sigma0 must be a square matrix of dimension length(beta0)')
  }

  ## assemble stan data
  standat = list(
    'N'           = nrow(X),
    'N0'          = nrow(X0),
    'p'           = ncol(X),
    'y'           = y,
    'y0'          = y0,
    'X'           = X,
    'X0'          = X0,
    'beta0'       = beta0,
    'Sigma0'      = Sigma0,
    'link'        = mu_link,
    'incl_offset' = incl_offset,
    'offset'      = offset,
    'offset0'     = offset0,
    'a0'          = a0
  )

  if ( family$family %in% c("gaussian", "Gamma", 'inverse.gaussian') ) {
    standat = c(standat, 'disp_shape' = disp.shape, 'disp_scale' = disp.scale)
  }

  ## call stan and return stanobject
  if (family$family == 'binomial') {
    return(
      rstan::sampling(
        object = stanmodels$glm_pp_fixed_bernoulli_post,
        data   = standat,
        ...
      )
    )
  }
  if ( family$family == "poisson" ) {
    return(
      rstan::sampling(
        object = stanmodels$glm_pp_fixed_poisson_post,
        data   = standat,
        ...
      )
    )
  }
  if ( family$family == "gaussian" ) {
    return(
      rstan::sampling(
        object = stanmodels$lm_pp_fixed_gaussian_post,
        data   = standat,
        ...
      )
    )
  }
  if ( family$family == "Gamma" ) {
    return(
      rstan::sampling(
        object = stanmodels$glm_pp_fixed_gamma_post,
        data   = standat,
        ...
      )
    )
  }
  if ( family$family == "inverse.gaussian" ) {
    return(
      rstan::sampling(
        object = stanmodels$glm_pp_fixed_invgaussian_post,
        data   = standat,
        ...
      )
    )
  }
  stop("Invalid family")
  return(NA)   ## never reached
}
