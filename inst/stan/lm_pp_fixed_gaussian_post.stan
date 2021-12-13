
data {
  int<lower=0>    N;              // number of obs current data
  int<lower=0>    N0;             // number of obs historical data
  int<lower=0>    p;              // number of regression coefficients
  vector[N]       y;              // response variable for current data
  vector[N0]      y0;             // response for historical data
  matrix[N,p]     X;              // des mat for current data
  matrix[N0,p]    X0;             // des mat for historical data
  vector[p]       beta0;          // initital prior mean for beta
  cov_matrix[p]   Sigma0;         // initial prior covariance for beta
  real<lower=0>   disp_shape;     // shape parameter for inv gamma prior for sigmasq
  real<lower=0>   disp_scale;     // scale parameter for inv gamma prior for sigmasq
  real<lower=0,upper=1>   a0;     // power prior parameter
  vector[N]       offset;         // offset for current data
  vector[N0]      offset0;        // offset for historical data
}
transformed data{
  vector[N] y_offset = y - offset;
  vector[N0] y0_offset = y0 - offset0;
}
parameters {
  vector[p] beta;
  real<lower=0> sigmasq;
}
model {
  real sigma = sqrt(sigmasq);
  beta    ~  multi_normal(beta0, sigmasq * Sigma0);              // conjugate init prior for beta
  sigmasq ~  inv_gamma(disp_shape, disp_scale);            // conjugate init prior for sigmasq
  target  += a0 * normal_id_glm_lpdf(y0_offset | X0, 0.0, beta, sigma); // power prior
  y_offset  ~ normal_id_glm(X, 0.0, beta, sigma);                       // likelihood
}

