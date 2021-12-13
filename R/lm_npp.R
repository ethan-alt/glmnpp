


#' Normalized power prior for for the normal linear model (identity link) using Stan
#'
#' Sample from the posterior of a normalized power prior (NPP) of a linear model (LM), where the normalizing constant is known.
#'
#' @param formula         an object of class \code{\link[stats]{formula}}.
#' @param data            a \code{\link[base]{data.frame}} of current data giving all variables in \code{formula}
#' @param histdata        a \code{\link[base]{data.frame}} of historical data giving all variables in \code{formula}
#' @param beta0           mean for initial prior on regression coefficients. Defaults to vector of 0s
#' @param Sigma0          covariance matrix for initial prior on regression coefficients. Defaults to \code{diag(100, ncol(X))}
#' @param offset          offset for current data. If \code{NULL}, no offset is utilized
#' @param offset0         offset for historical data. If \code{NULL}, no offset is utilized
#' @param sigmasq.shape   shape parameter for inverse-gamma prior on variance parameter. Defaults to \code{2.1}
#' @param sigmasq.scale   scale parameter for inverse-gamma prior on variance parameter. Defaults to \code{1.1}
#' @param prec.shape      shape parameter for gamma prior on precision parameter (inverse of \code{sigmasq}). Ignored if sigmasq.shape is specified
#' @param prec.rate       rate parameter for gamma prior on precision parameter (inverse of \code{sigmasq}). Ignored if sigmasq.scale is specified
#' @param a0.shape1       shape1 parameter for beta prior on the power prior parameter. When `a0.shape1==1` and `a0.shape2==1`, a uniform prior is utilized
#' @param a0.shape2       shape2 parameter for beta prior on the power prior parameter. When `a0.shape1==1` and `a0.shape2==1`, a uniform prior is utilized
#' @param ...             optional parameters to pass onto `rstan::sampling`
#'
#' @return an object of class [rstan::stanfit] returned by `rstan::sampling`
#' @examples
#' ## Generate current and historical data
#' set.seed(123)
#' N = 100
#' N0 = 50
#' x = rnorm(N)
#' y = 1 + 0.5 * x + rnorm(N)
#' x0 = rnorm(N0)
#' y0 = 1 + 0.5 * x0 + rnorm(N0)
#' data = data.frame('y' = y, 'x' = x)
#' histdata = data.frame('y' = y0, 'x' = x0)
#'
#' ## Perform MCMC of normalized power prior for normal linear model
#' fit = lm_npp(y ~ x, data, histdata)
#'
#' @export
lm_npp = function(
  formula, data, histdata, beta0 = NULL, Sigma0 = NULL, offset = NULL,
  offset0 = NULL, sigmasq.shape = NULL, sigmasq.scale = NULL,
  prec.shape = NULL, prec.rate = NULL, a0.shape1 = 1, a0.shape2 = 1, ...
) {
  ## make sure formula is two-sided
  if(!(formula.tools::is.two.sided(formula))) {
    stop("formula must be two-sided")
  }
  if(any(is.na(data)))     { stop('There are missing values in data') }
  if(any(is.na(histdata))) { stop('There are missing values in histdata') }

  ## Check hyperparameters for sigma^2
  if (!is.null(prec.rate) & is.null(sigmasq.scale) )
    sigmasq.scale = prec.rate
  if (!is.null(prec.shape) & is.null(sigmasq.shape) )
    sigmasq.shape = prec.shape
  if ( is.null(sigmasq.shape) )
    sigmasq.shape = 2.1
  if(is.null(sigmasq.scale))
    sigmasq.scale = 1.1
  ## get design matrices
  X  = model.matrix(formula, data)
  X0 = model.matrix(formula, histdata)

  ## get responses
  y  = data[, all.vars(formula)[1]]
  y0 = histdata[, all.vars(formula)[1]]

  ## check if offset is included
  if ( !is.null(offset) ) {
    ## check dimensions of offset
    if ( length(offset) != length(y) )   { stop("length of offset must match data") }
    if ( length(offset0) != length(y0) ) { stop("length of offset0 must match histdata") }
  } else {
      offset = rep(0, length(y))
      offset0 = rep(0, length(y0))
  }

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
    'N'             = nrow(X),
    'N0'            = nrow(X0),
    'p'             = ncol(X),
    'y'             = y,
    'y0'            = y0,
    'X'             = X,
    'X0'            = X0,
    'mu0'           = beta0,
    'Sigma0'        = Sigma0,
    'sigmasq_shape' = sigmasq.shape,
    'sigmasq_scale' = sigmasq.scale,
    'a0_shape1'     = a0.shape1,
    'a0_shape2'     = a0.shape2,
    'offset'        = offset,
    'offset0'       = offset0
  )

  ## call stan and return stanobject
  rstan::sampling(
    object = stanmodels$lm_npp_gaussian,
    data   = standat,
    ...
  )
}
