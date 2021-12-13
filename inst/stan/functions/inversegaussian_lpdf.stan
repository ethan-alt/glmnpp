real inversegaussian_lpdf(vector y, vector mu, real lambda) {
  return 0.5 * (
      rows(y) * log(lambda / (6.283185307179586)) - 3 * sum(log(y))
    - lambda * dot_self( (y - mu) .* inv(mu .* sqrt(y)) )
  );
}
