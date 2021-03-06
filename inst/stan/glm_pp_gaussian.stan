
functions {
#include /functions/link.stan
}


data {
  int<lower=0>          nobs;          // number of observations
  int<lower=0>          p;             // number of predictors
  real                  y0[nobs];      // historical data
  matrix[nobs,p]        X;             // design matrix
  real<lower=0,upper=1> a0;            // power prior parameter
  vector[p]             beta0;         // initial prior mean of beta
  matrix[p,p]           Sigma0;        // initial prior covariance matrix of beta
  real<lower=0>         disp_shape;    // shape parameter for inverse-gamma prior on dispersion
  real<lower=0>         disp_scale;    // scale parameter for inverse-gamma prior on dispersion
  int<lower=1,upper=9>  link;          // index of link function
  int<lower=0,upper=1>  incl_offset;   // whether an offset is included
  vector[nobs]          offset;        // offset (defaults to vector of 0s in R)
}
transformed data {
  vector[nobs] y0vec = to_vector(y0);
  y0vec = y0vec - offset;
}

// p+1 params: p-dim vector of regression coefficients and scalar inverse dispersion
parameters {
  vector[p] beta;
  real<lower=0> dispersion;
}

// Assume beta is a priori MVN; obtain posterior
// based on power prior
model {
  vector[( (link != 1) || (incl_offset == 1) ) ? nobs :  0] eta;
  vector[(link != 1) ? nobs :  0] mu;
  real sigma = sqrt(dispersion);

  // initial priors
  beta       ~ multi_normal(beta0, Sigma0);        // MVN initial prior on beta
  dispersion ~ inv_gamma(disp_shape, disp_scale);  // inverse-gamma initial prior on dispersion

  if ( a0 > 0 ) {
    if ( link == 1 )
      target += a0 * normal_id_glm_lpdf(y0vec | X, 0.0, beta, sigma);
    else {
      eta = X * beta;
      mu      = linkinv(eta, link);
      target += a0 * normal_lpdf(y0vec | mu, sigma);
    }
  }
}
