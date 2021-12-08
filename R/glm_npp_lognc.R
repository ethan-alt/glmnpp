


#' Estimate logarithm of the normalizing constant for the power prior
#'
#' Estimate the logarithm of the normalizing constant via importance sampling or bridge sampling.
#' For importance sampling, the function uses the \code{optimizing} function in Stan to compute the prior mode and Hessian, then uses a normal approximation.
#' Otherwise, MCMC sampling is performed and bridge sampling is utilized to estimate the logarithm of the normalizing constant.
#' When a dispersion parameter is unknown, it is better to use bridge sampling (\code{method = 'bridge'}).
#'
#' @include glm_npp_prior.R
#'
#' @param formula         an object of class \code{\link[stats]{formula}}.
#' @param family          an object of type \code{\link[stats]{family}} giving distribution and link function
#' @param histdata        a \code{\link[base]{data.frame}} of historical data giving all variables in \code{formula}
#' @param a0              a vector of values between 0 and 1 giving values of the normalizing constant to compute. One of \code{a0} or \code{a0.n} must be specified
#' @param a0.n            a positive integer giving the number of (equally spaced) values to compute the normalizing constant of the power prior. One of \code{a0} or \code{a0.n} must be specified
#' @param beta0           mean for initial prior on regression coefficients. Defaults to vector of 0s
#' @param Sigma0          covariance matrix for initial prior on regression coefficients. Defaults to \code{diag(100, ncol(X))}
#' @param offset          offset in GLM. If \code{NULL}, no offset is utilized
#' @param disp.shape      shape parameter for inverse-gamma prior on dispersion parameter (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param disp.scale      rate parameter for inverse-gamma prior on dispersion parameter (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param method          character vector giving which method to use for importance sampling. Acceptable values are `"bridge"` or`"importance"`, corresponding to bridge sampling and importance sampling, respectively
#' @param nsmpl           (optional) number of importance samples to take (ignored if \code{method == 'bridge'})
#' @param bridge.args     (optional) parameters to pass onto `bridgesampling::bridge_sampler` (otherwise, default is performed) (ignored if \code{method == 'importance'})
#' @param ...             (optional) parameters to pass onto `rstan::sampling` (ignored if \code{method == 'importance'})
#'
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
#' ## Estimate logarithm of normalizing constant using bridge sampling
#' bridge = glm_npp_lognc(y ~ x, binomial(), histdata, a0 = c(0.10, 0.50), method = 'bridge')
#'
#' ## Estimate logarithm of normalizing constant using importance sampling
#' ##   Note: for demonstration purposes only. nsmpl should be at least 10000 in practice
#' importance = glm_npp_lognc(y ~ x, binomial(), histdata, a0 = c(0.10, 0.50), method = 'importance', nsmpl = 2000)
#'
#' ## Compare results
#' cbind('a0' = bridge[, 1], 'bridge' = bridge[, 2], 'importance' = importance[, 2])
#'
#'
#' @return a \code{data.frame} giving a0 and the logarithm of the normalizing constant evaluted at a0
#'
#' @export
glm_npp_lognc = function(formula, family, histdata, a0 = NULL, a0.n = NULL, beta0 = NULL, Sigma0 = NULL, offset = NULL, disp.shape = 2.1, disp.scale = 1.1, method = 'bridge', nsmpl = 1000, bridge.args = NULL, ...) {
  ## Check if a0 or a0.n are specified
  if( !(is.null(a0)) & !(is.null(a0.n)) ) {
    warning("a0 and a0.n both specified. Using the input for a0.")
    a0.n = NULL
  }
  if ( is.null(a0) ) {
    if(is.null(a0.n)) stop("One of a0.n or a0 must be specified")
    a0   = seq(0, 1, length.out = a0.n)
    a0.n = NULL
  }

  ## Supply warning if linear model with identity link
  if ( family$family == 'gaussian' & family$link == 'identity' )
    message("Estimating the normalizing constant is not necessary for gaussian models with identity link functions. Use the function lm_npp for better accuracy.")

  ## check method
  if (!(method %in% c('bridge', 'importance'))) stop('method must be "bridge" or "importance"')

  ## Initialize result
  lognc = numeric(length(a0))
  if ( method == 'bridge' ) {
    for ( i in seq_along(a0) ) {
      lognc[i] = glm_npp_lognc_bridge(formula, family, histdata, a0[i], beta0, Sigma0, offset, disp.shape, disp.scale, bridge.args, ...)
    }
  }
  if ( method == 'importance' ){
    for ( i in seq_along(a0) ) {
      lognc[i] = glm_npp_lognc_importance(formula, family, histdata, a0[i], beta0, Sigma0, offset, disp.shape, disp.scale, nsmpl = nsmpl)
    }
  }
  ## result is a data.frame giving (a0, lognc)
  a0lognc                 = data.frame('a0' = a0, 'lognc' = lognc)
  a0lognc                 = a0lognc[order(a0lognc$a0), ]
  attr(a0lognc, 'method') = method
  return(a0lognc)
}




#' Estimate logarithm of normalizing constant using bridge sampling
#'
#' Estimate the logarithm of the normalizing constant via importance sampling. Uses the \code{optimizing} function in Stan to compute
#' the posterior mode and Hessian, then uses a normal approximation
#'
#' @param formula         an object of class \code{\link[stats]{formula}}.
#' @param family          an object of type \code{\link[stats]{family}} giving distribution and link function
#' @param histdata        a \code{\link[base]{data.frame}} of historical data giving all variables in \code{formula}
#' @param a0              a positive scalar no larger than 1 giving the power prior parameter
#' @param beta0           mean for initial prior on regression coefficients. Defaults to vector of 0s
#' @param Sigma0          covariance matrix for initial prior on regression coefficients. Defaults to \code{diag(100, ncol(X))}
#' @param offset          offset in GLM. If \code{NULL}, no offset is utilized
#' @param disp.shape      shape parameter for inverse-gamma prior on dispersion parameter (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param disp.scale      rate parameter for inverse-gamma prior on dispersion parameter (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param bridge.args     (optional) parameters to pass onto `bridgesampling::bridge_sampler` (otherwise, default is performed)
#' @param ...             (optional) parameters to pass onto `rstan::sampling`
glm_npp_lognc_bridge = function(formula, family, histdata, a0, beta0, Sigma0, offset, disp.shape, disp.scale, bridge.args = NULL, ...) {
  fit = glm_npp_prior(
    formula, family, histdata = histdata, a0 = a0, beta0 = beta0, Sigma0 = Sigma0, offset = offset, disp.shape = disp.shape, disp.scale = disp.scale, ...
  )
  ## Check if effective sample size / rhat is OK
  diag = rstan::summary(fit)$summary[, c('n_eff', 'Rhat')]
  rhat = max(diag[, 2])
  neff = min(diag[, 1])
  if ( rhat > 1.05 )
    stop("Maximum rhat is larger than 1.05. Try increasing number of samples")
  if ( neff < 1000 )
    stop("Minimum effective sample size is less than 1000. Try increasing number of samples")
  if ( is.null(bridge.args) ) {
    res = bridgesampling::bridge_sampler(fit)
  } else {
    args = c('samples' = fit, bridge.args)
    res = do.call(bridgesampling::bridge_sampler, args)
  }
  return(res$logml)
}






#' Estimate logarithm of normalizing constant using importance sampling
#'
#' Estimate the logarithm of the normalizing constant via importance sampling. Uses the \code{optimizing} function in Stan to compute
#' the posterior mode and Hessian, then uses a normal approximation
#'
#' @param formula         an object of class \code{\link[stats]{formula}}.
#' @param family          an object of type \code{\link[stats]{family}} giving distribution and link function
#' @param histdata        a \code{\link[base]{data.frame}} of historical data giving all variables in \code{formula}
#' @param a0              a positive scalar no larger than 1 giving the power prior parameter
#' @param beta0           mean for initial prior on regression coefficients. Defaults to vector of 0s
#' @param Sigma0          covariance matrix for initial prior on regression coefficients. Defaults to \code{diag(100, ncol(X))}
#' @param offset          offset in GLM. If \code{NULL}, no offset is utilized
#' @param disp.shape      shape parameter for inverse-gamma prior on dispersion parameter (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param disp.scale      rate parameter for inverse-gamma prior on dispersion parameter (for Gaussian and gamma models). Ignored for binomial and Poisson models
#' @param nsmpl           number of importance samples to take
#'
#' @return scalar giving log normalizing constant
glm_npp_lognc_importance = function(formula, family, histdata, a0, beta0, Sigma0, offset, disp.shape, disp.scale, nsmpl) {
  ## get design matrix
  X = model.matrix(formula, histdata)

  ## get response
  y0 = histdata[, all.vars(formula)[1]]


  ## get mu link function as integer
  links = c('identity', 'log', 'logit', 'inverse', 'probit', 'cauchit', 'cloglog', 'sqrt', '1/mu^2')
  mu_link = which(links == family$link)[1]
  if ( length(mu_link) == 0 ) { stop(paste('Link must be one of', paste(links, collapse = ', '))) }

  ## get offset if applicable
  incl_offset = 1
  if ( is.null(offset) ) {
    incl_offset = 0
    offset = rep(0, length(y0))
  }

  ## default mean for beta is 0
  if(is.null(beta0)) {
    beta0 = rep(0, ncol(X))
  }

  ## default covariance for beta is diag(100)
  if(is.null(Sigma0)) {
    Sigma0 = diag(100, ncol(X))
  }


  ## assemble stan data
  standat = list(
    'nobs'        = nrow(X),
    'p'           = ncol(X),
    'y0'          = y0,
    'X'           = X,
    'a0'          = a0,
    'beta0'       = beta0,
    'Sigma0'      = Sigma0,
    'link'        = mu_link,
    'incl_offset' = incl_offset,
    'offset'      = offset
  )

  if ( family$family %in% c("gaussian", "Gamma") ) {
    standat = c(standat, 'disp_shape' = disp.shape, 'disp_scale' = disp.scale)
  }

  ## call stan and return stanobject
  if (family$family == 'binomial') {
    opt = rstan::optimizing(
      object  = stanmodels$glm_pp_bernoulli,
      data    = standat,
      hessian = T
    )
    stanobj = suppressMessages(
      rstan::sampling(
        object = stanmodels$glm_pp_bernoulli,
        data   = standat,
        chains = 0
      )
    )
  }

  if ( family$family == "poisson" ) {
    opt = rstan::optimizing(
      object = stanmodels$glm_pp_poisson,
      data    = standat,
      hessian = T
    )
    stanobj = suppressMessages(
      rstan::sampling(
        object = stanmodels$glm_pp_poisson,
        data   = standat,
        chains = 0
      )
    )
  }

  if ( family$family == "gaussian" ) {
    opt = rstan::optimizing(
      object  = stanmodels$glm_pp_gaussian,
      data    = standat,
      hessian = T
    )
    stanobj = suppressMessages(
      rstan::sampling(
        object  = stanmodels$glm_pp_gaussian,
        data   = standat,
        chains = 0
      )
    )
  }

  if ( family$family == "Gamma" ) {
    opt = rstan::optimizing(
      object  = stanmodels$glm_pp_gamma,
      data    = standat,
      hessian = T
    )
    stanobj = suppressMessages(
      rstan::sampling(
        object  = stanmodels$glm_pp_gamma,
        data   = standat,
        chains = 0
      )
    )
  }

  ## obtain importance sample from MVN
  invneghess = chol2inv(chol(-opt$hessian) )
  smpl = mvtnorm::rmvnorm(n = nsmpl, mean = opt$par, sigma = invneghess )

  ## compute log weights
  logW = apply(smpl, 1, function(x) rstan::log_prob(stanobj, x) )
  logW = logW - mvtnorm::dmvnorm(smpl, opt$par, invneghess, log = T)

  ## log-sum-exp trick to compute log NC
  M     = max(logW)
  lognc = -log(nsmpl) + M + log( sum( exp( logW - M ) ) )

  return(lognc)
}
