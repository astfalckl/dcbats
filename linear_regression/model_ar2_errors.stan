data {
  int<lower=1> T;
  int<lower=1> p;
  matrix[T, p] Z;
  vector[T] y;
  real<lower=0> temper;
}

parameters {
  real alpha;
  vector[p] beta;
  real phi1;
  real phi2;
  real<lower=0> sigma2;
  real eps_m1;   // epsilon_{-1}
  real eps_0;    // epsilon_0
}

transformed parameters {
  real<lower=0> sigma;
  vector[T] r;

  sigma = sqrt(sigma2);
  r = y - alpha - Z * beta;
}

model {
  // priors
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);
  phi1  ~ normal(0, 10);
  phi2  ~ normal(0, 10);
  sigma2 ~ inv_gamma(3, 10);

  eps_m1 ~ normal(0, sigma);
  eps_0  ~ normal(0, sigma);

  // tempered likelihood
  target += temper * normal_lpdf(r[1] | phi1 * eps_0 + phi2 * eps_m1, sigma);

  if (T >= 2) {
    target += temper * normal_lpdf(r[2] | phi1 * r[1] + phi2 * eps_0, sigma);
  }

  for (t in 3:T) {
    target += temper * normal_lpdf(r[t] | phi1 * r[t - 1] + phi2 * r[t - 2], sigma);
  }
}

