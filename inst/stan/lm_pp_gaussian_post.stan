
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
transformed data{
  matrix[K,K] hist_prec       = crossprod(X0); // X0'X0
  vector[K]   hist_mle        = mdivide_left_spd(hist_prec, X0' * y0); // betahat0
  matrix[K,K] prior_prec      = inverse_spd(Sigma0);   // prior precision matrix
  vector[K]   hist_prec_mle   = hist_prec * hist_mle;  // X0'X0 * betahat0
  vector[K]   prior_prec_mean = prior_prec * mu0;      // priorPrec * mu0
  real        sumy0sq         = dot_self(y0);
  real        quad_b0         = quad_form(prior_prec, mu0);  // mu0' Omega0 mu0
}
parameters {
  vector[K] beta;
  real<lower=0> sigmasq;
  real<lower=0,upper=1> a0;
}
model {
  real sigma = sqrt(sigmasq);
  // compute prior mean and precision given a0
  matrix[K,K] prior_prec_a0 = a0 * hist_prec + prior_prec;
  vector[K]   prior_mean_a0 =
    mdivide_left_spd(prior_prec_a0, a0 * hist_prec_mle + prior_prec_mean );

  // compute prior parameters for sigmasq given a0
  real sigmasq_shape_a0 = sigmasq_shape + 0.5 * a0 * N0;
  real sigmasq_rate_a0 = sigmasq_rate +
    0.5 * ( a0 * sumy0sq + quad_b0 - quad_form(prior_prec_a0, prior_mean_a0) );

  // prior given a0
  beta       ~ multi_normal_prec(prior_mean_a0, inv(sigmasq) * prior_prec_a0);
  sigmasq    ~ inv_gamma(sigmasq_shape_a0, sigmasq_rate_a0);
  if ( a0_shape1 != 1 || a0_shape2 != 1 )
    a0 ~ beta(a0_shape1, a0_shape2);

  // likelihood of current data and power prior
  y ~ normal_id_glm(X, 0.0, beta, sigma);
}

