#' Normalized Power Prior for Generalized Linear Models
#'
#' @description User-friendly package to implement the normalized power prior
#'              for generalized linear models. Contains a function that samples
#'              from the posterior distribution of a normalized power prior
#'              directly for the normal linear model. For generalized linear
#'              models, we use a two-step process that first estimates the
#'              normalizing constant over a fine grid of values, then uses
#'              those estimates to sample from the posterior distribution.
#'
#' @docType package
#' @name glmnpp-package
#' @aliases glmnpp
#' @useDynLib glmnpp, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import stats
#' @importFrom rstan sampling
#' @importFrom formula.tools is.two.sided
#' @importFrom mvtnorm rmvnorm
#' @importFrom bridgesampling bridge_sampler
#' @importFrom utils capture.output
#'
#' @references
#' Carvallho, L. M. and Ibrahim, J. G. (2021): On the normalized power prior. Statistics in Medicine. doi:10.1002/sim.9124
#'
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
NULL
