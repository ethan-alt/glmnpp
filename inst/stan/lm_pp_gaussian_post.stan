
data {
  int<lower=0>    N;             // number of obs current data
  int<lower=0>    N0;            // number of obs historical data
  int<lower=0>    K;             // number of coeffs incl. intercept
  vector[N]       y;             // response variable for current data
  vector[N0]      y0;            // response for historical data
  matrix[N,K]     X;             // des mat for current data
  matrix[N0,K]    X0;            // des mat for historical data
  vector[K]       mu0;           // initital prior mean for beta
  cov_matrix[K]   Sigma0;        // initial prior covariance for beta
  real<lower=0>   sigmasq_shape; // shape parameter for sigmasq
  real<lower=0>   sigmasq_rate;  // rate parameter for sigmasq
  real<lower=0>   a0_shape1;     // shape parameter for a0
  real<lower=0>   a0_shape2;     // shape2 parameter for a0
}
parameters {
  vector[K] beta;
  real<lower=0> sigmasq;
  real<lower=0,upper=1> a0;
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // (initial) priors for parameters
  real sigma = sqrt(sigmasq);
  beta       ~ multi_normal(mu0, Sigma0);
  sigmasq    ~ inv_gamma(sigmasq_shape, sigmasq_rate);
  if ( a0_shape1 != 1 || a0_shape2 != 1 )
    a0 ~ beta(a0_shape1, a0_shape2);
  
  // likelihood of current data and power prior  
  y ~ normal_id_glm_lpdf(X, 0.0, beta, sigma);
  target += normal_id_glm_lpdf(y0 | X0, 0.0, beta, sigma * inv_sqrt(a0));
}

