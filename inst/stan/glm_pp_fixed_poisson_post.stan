
functions {
#include /functions/link.stan
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0>          N;             // number of observations in current data
  int<lower=0>          N0;            // number of observations in historical data
  int<lower=0>          p;             // number of predictors
  int<lower=0>          y[N];          // current response vector
  int<lower=0>          y0[N0];        // historical response vector
  matrix[N,p]           X;             // design matrix for current data
  matrix[N0,p]          X0;            // design matrix for historical data
  vector[p]             beta0;         // initial prior mean of beta
  matrix[p,p]           Sigma0;        // initial prior covariance matrix of beta
  int<lower=1,upper=9>  link;          // index of link function
  int<lower=0,upper=1>  incl_offset;   // whether an offset is included
  vector[N]             offset;        // offset for current data (defaults to vector of 0s in R)
  vector[N0]            offset0;       // offset for historical data (defaults to vector of 0s in R)
  real<lower=0,upper=1> a0;            // power prior parameter
}

// Two parameters: regression coefficient and power prior param
parameters {
  vector[p] beta;
}

model {
  // obtain linear predictor (adding offset if applicable)
  vector[( (link != 2) || (incl_offset == 1) ) ? N :  0] eta;
  vector[( (link != 2) || (incl_offset == 1) ) ? N0 :  0] eta0;
  vector[(link != 2) ? N :  0] mu;
  vector[(link != 2) ? N0 :  0] mu0;

  // add to target initial prior
  beta    ~ multi_normal(beta0, Sigma0);

  // add log likelihood and power prior--using the best stan function available per scenario

  // if canonical link used and no offset
  if ( link == 2 && incl_offset == 0 ) {
    target += a0 * poisson_log_glm_lpmf(y0 | X0, 0, beta);
    target += poisson_log_glm_lpmf(y | X, 0, beta);
  }
  else {
    eta  = X  * beta;
    eta0 = X0 * beta0;
    if ( incl_offset == 1 ){
      eta  += offset;
      eta0 += offset0;
    }
    if ( link == 2 ) {    // if canonical link used with offset
      target += a0 * poisson_log_lpmf(y0 | eta0);
      target += poisson_log_lpmf(y | eta);
    }
    else {                 // non-canonical link used, possibly with offset
      mu  = linkinv(eta, link);
      mu0 = linkinv(eta0, link);
      target += a0 * poisson_lpmf(y0 | mu0);
      target += poisson_lpmf(y | mu);
    }
  }
}
