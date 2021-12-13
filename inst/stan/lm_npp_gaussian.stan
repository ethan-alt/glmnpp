
data {
  int<lower=0>    N;              // number of obs current data
  int<lower=0>    N0;             // number of obs historical data
  int<lower=0>    p;              // number of coeffs incl. intercept
  vector[N]       y;              // response variable for current data
  vector[N0]      y0;             // response for historical data
  matrix[N,p]     X;              // des mat for current data
  matrix[N0,p]    X0;             // des mat for historical data
  vector[p]       mu0;            // initital prior mean for beta
  cov_matrix[p]   Sigma0;         // initial prior covariance for beta
  real<lower=0>   sigmasq_shape;  // shape parameter for inv gamma prior for sigmasq
  real<lower=0>   sigmasq_scale;  // scale parameter for inv gamma prior for sigmasq
  real<lower=0>   a0_shape1;      // shape parameter for a0
  real<lower=0>   a0_shape2;      // shape2 parameter for a0
  vector[N]       offset;         // offset for current data
  vector[N0]      offset0;        // offset for historical data
}
transformed data{
  vector[N]   y_offset        = y - offset;
  vector[N0]  y0_offset       = y0 - offset0;
  matrix[p,p] hist_prec       = crossprod(X0); // X0'X0
  vector[p]   hist_mle        = mdivide_left_spd(crossprod(X0), X0' * y0_offset); // betahat0
  matrix[p,p] prior_prec      = inverse_spd(Sigma0);   // prior precision matrix
  vector[p]   hist_prec_mle   = hist_prec * hist_mle;  // X0'X0 * betahat0
  vector[p]   prior_prec_mean = prior_prec * mu0;      // Omega0 * beta0
  real        sumy0sq         = dot_self(y0_offset);          // sum(y0^2)
  real        quad_b0         = quad_form(prior_prec, mu0);  // mu0' Omega0 mu0
}
parameters {
  vector[p] beta;
  real<lower=0> sigmasq;
  real<lower=0,upper=1> a0;
}
model {
  // Compute parameters for implied normal-inverse-Gamma prior on beta, sigma^2
    // Omega* =  a0 X0'X0 + Omega0 = precision matrix | sigmasq, a0
    matrix[p,p] prec_a0  = a0 * hist_prec + prior_prec;
    // mu* = (a0 X0'X0 + Omega0)^(-1) a_0 X0'X0 betahat0 + Omega0 beta0 = mean | sigmasq, a0
    vector[p] mean_a0  = mdivide_left_spd(prec_a0, a0 * hist_prec_mle + prior_prec_mean);
    // shape for sigmasq = alpha0 + 0.5 * a0 * N0
    real shape_a0 = sigmasq_shape + 0.5 * a0 * N0;
    // scale for sigmasq = eta0 + a0 * y0' y0 + beta0' Omega0 beta0 - mu*' Omega* mu*
    real scale_a0 = sigmasq_scale + 0.5 *
      ( a0 * sumy0sq + quad_b0 - quad_form(prec_a0, mean_a0) );

  // likelihood of current data
  y_offset ~ normal_id_glm(X, 0.0, beta, sqrt(sigmasq));

  // power prior is normal-inverse Gamma
  target += inv_gamma_lpdf(sigmasq | shape_a0, scale_a0);
  target += multi_normal_prec_lpdf(beta | mean_a0, prec_a0 * inv(sigmasq) );
}

