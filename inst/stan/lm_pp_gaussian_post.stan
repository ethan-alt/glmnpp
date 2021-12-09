
data {
  int<lower=0>    N;              // number of obs current data
  int<lower=0>    N0;             // number of obs historical data
  int<lower=0>    K;              // number of coeffs incl. intercept
  vector[N]       y;              // response variable for current data
  vector[N0]      y0;             // response for historical data
  matrix[N,K]     X;              // des mat for current data
  matrix[N0,K]    X0;             // des mat for historical data
  vector[K]       mu0;            // initital prior mean for beta
  cov_matrix[K]   Sigma0;         // initial prior covariance for beta
  real<lower=0>   sigmasq_shape;  // shape parameter for inv gamma prior for sigmasq
  real<lower=0>   sigmasq_scale;  // scale parameter for inv gamma prior for sigmasq
  real<lower=0>   a0_shape1;      // shape parameter for a0
  real<lower=0>   a0_shape2;      // shape2 parameter for a0
}
transformed data{
  matrix[K,K] hist_prec       = crossprod(X0); // X0'X0
  vector[K]   hist_mle        = mdivide_left_spd(crossprod(X0), X0' * y0); // betahat0
  matrix[K,K] prior_prec      = inverse_spd(Sigma0);   // prior precision matrix
  vector[K]   hist_prec_mle   = hist_prec * hist_mle;  // X0'X0 * betahat0
  vector[K]   prior_prec_mean = prior_prec * mu0;      // priorPrec * mu0
  real        sumy0sq         = dot_self(y0);          // sum(y0^2)
  real        quad_b0         = quad_form(prior_prec, mu0);  // mu0' Omega0 mu0
}
parameters {
  vector[K] beta;
  real<lower=0> sigmasq;
  real<lower=0,upper=1> a0;
}
model {
  real sigma = sqrt(sigmasq);
  matrix[K,K] prec_a0 = a0 * hist_prec + prior_prec;  // prior prec using npp
  vector[K]   mean_a0 = mdivide_left_spd(prec_a0, a0 * hist_prec_mle + prior_prec_mean);
    // prior mean using npp

  real sigmasq_shape_a0 = sigmasq_shape + 0.5 * (a0 * N0);          // a0 + 0.5 * (a0 * n0 - p)
  real sigmasq_scale_a0 = sigmasq_scale +                           // b0 + 0.5 * (a0 * y0'y0 + beta0'Omega0 beta0 - mu*' Omega* mu*)
    0.5 * (a0 * sumy0sq + quad_b0 - quad_form(prec_a0, mean_a0));

  // likelihood of current data and power prior
  y ~ normal_id_glm(X, 0.0, beta, sigma);
  target += a0 * normal_id_glm_lpdf(y0 | X0, 0.0, beta, sigma);

  // priors
  sigmasq ~ inv_gamma(sigmasq_shape, sigmasq_scale);
  beta    ~ multi_normal_prec(mu0, prior_prec * inv(sigmasq));

  // add normalizing constant
  target += -0.5 * log_determinant(a0 * hist_prec + prior_prec); // log | a0 X0'X0 + Omega0 |
  target += -lgamma(sigmasq_shape_a0);
  target += sigmasq_shape_a0 * log(sigmasq_scale_a0);
}

