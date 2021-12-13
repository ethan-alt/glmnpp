
functions {
#include /functions/link.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0>          N;             // number of observations in current data
  int<lower=0>          N0;            // number of observations in historical data
  int<lower=0>          p;             // number of predictors
  real<lower=0>         y[N];          // current response vector
  real<lower=0>         y0[N0];        // historical response vector
  matrix[N,p]           X;             // design matrix for current data
  matrix[N0,p]          X0;            // design matrix for historical data
  vector[p]             beta0;         // initial prior mean of beta
  matrix[p,p]           Sigma0;        // initial prior covariance matrix of beta
  int<lower=1,upper=9>  link;          // index of link function
  int<lower=0,upper=1>  incl_offset;   // whether an offset is included
  vector[N]             offset;        // offset for current data (defaults to vector of 0s in R)
  vector[N0]            offset0;       // offset for historical data (defaults to vector of 0s in R)
  real<lower=0,upper=1> a0;            // power prior parameter
  real<lower=0>         disp_shape;    // shape parameter for inverse-gamma prior on inverse dispersion
  real<lower=0>         disp_scale;    // rate parameter on inverse-gamma prior for inverse dispersion
}

// Two parameters: regression coefficient and power prior param
parameters {
  vector[p] beta;
  real<lower=0> dispersion;
}

model {
  real alpha = inv(dispersion);
  vector[N]  mu;
  vector[N0] mu0;
  // obtain linear predictor (adding offset if applicable)
  vector[N0] eta0 = X0 * beta;
  vector[N]  eta  = X  * beta;
  if (incl_offset == 1) {
    eta0 = eta0 + offset0;
    eta  = eta  + offset;
  }

  // compute means
  mu  = linkinv(X * beta, link);
  mu0 = linkinv(X0 * beta, link);

  // prior for beta is MVN and for dispersion is inverse-gamma


  beta       ~ multi_normal(beta0, Sigma0);                      // initial prior on beta is MVN
  dispersion ~ inv_gamma(disp_shape, disp_scale);                // initial prior on dispersion is inverse-gamma
  target     += a0 * gamma_lpdf(y0 | alpha, alpha * inv(mu0));   // power prior
  target     += gamma_lpdf(y | alpha, alpha * inv(mu));          // likelihood
}
