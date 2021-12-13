


#' Posterior inference using a power prior with fixed a0
#' for the normal linear model
#'
#' Sample from the posterior of a power prior (PP) of a linear model (LM).
#' This is a wrapper for the function `glm_pp` with
#' `family == gaussian('identity')`.
#'
#' @param formula         an object of class \code{\link[stats]{formula}}.
#' @param data            a \code{\link[base]{data.frame}} of current data giving all variables in \code{formula}
#' @param histdata        a \code{\link[base]{data.frame}} of historical data giving all variables in \code{formula}
#' @param a0              power prior parameter (between 0 and 1)
#' @param beta0           mean for initial prior on regression coefficients. Defaults to vector of 0s
#' @param Sigma0          covariance matrix for initial prior on regression coefficients. Defaults to \code{diag(100, ncol(X))}
#' @param offset          offset for current data. If \code{NULL}, no offset is utilized
#' @param offset0         offset for historical data. If \code{NULL}, no offset is utilized
#' @param sigmasq.shape   shape parameter for inverse-gamma prior on variance parameter. Defaults to \code{2.1}
#' @param sigmasq.scale   scale parameter for inverse-gamma prior on variance parameter. Defaults to \code{1.1}
#' @param prec.shape      shape parameter for gamma prior on precision parameter (inverse of \code{sigmasq}). Ignored if sigmasq.shape is specified
#' @param prec.rate       rate parameter for gamma prior on precision parameter (inverse of \code{sigmasq}). Ignored if sigmasq.scale is specified
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
#' ## Perform MCMC of non-normalized power prior for normal linear model
#' fit = lm_pp(y ~ x, data, histdata, a0 = 0.5)
#'
#' @export
lm_pp = function(
  formula, data, histdata, a0 = 0.5, beta0 = NULL, Sigma0 = NULL, offset = NULL,
  offset0 = NULL, sigmasq.shape = NULL, sigmasq.scale = NULL,
  prec.shape = NULL, prec.rate = NULL, ...
) {
  ## Check hyperparameters for sigma^2
  if (!is.null(prec.rate) & is.null(sigmasq.scale) )
    sigmasq.scale = prec.rate
  if (!is.null(prec.shape) & is.null(sigmasq.shape) )
    sigmasq.shape = prec.shape
  if ( is.null(sigmasq.shape) )
    sigmasq.shape = 2.1
  if(is.null(sigmasq.scale))
    sigmasq.scale = 1.1

  return(
    glm_pp(
      formula, family = gaussian('identity'),
      data = data, histdata = histdata, a0 = a0, beta0 = beta0, Sigma0 = Sigma0,
      offset = offset, offset0 = offset0,
      disp.shape = sigmasq.shape, disp.scale = sigmasq.scale,
      ...
    )
  )
}
