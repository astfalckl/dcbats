data {
  int<lower=2> T;
  int<lower=2> dim;
  array[T] vector[dim] Y;
  real<lower=0> power;
  real<lower=0> sigma1_0;
  real<lower=0> sigma2_0;
}

transformed data {
  if (dim != 2) {
    reject("This model is hard-coded for dim = 2.");
  }
}

parameters {
  real<lower=0, upper=1> a1_raw;
  real<lower=0, upper=1> b1_raw;
  real<lower=0, upper=1> a2_raw;
  real<lower=0, upper=1> b2_raw;

  real<lower=0> w1;
  real<lower=0> w2;
  real<lower=0, upper=1> r;
  real mu1;
  real mu2;
}

transformed parameters {
  real<lower=0, upper=1> a1 = a1_raw;
  real<lower=0, upper=1 - a1> b1 = (1 - a1) * b1_raw;

  real<lower=0, upper=1> a2 = a2_raw;
  real<lower=0, upper=1 - a2> b2 = (1 - a2) * b2_raw;
}

model {
  vector[T] sigma1;
  vector[T] sigma2;
  matrix[2, 2] Sigma;
  vector[2] Mu;

  sigma1[1] = sigma1_0;
  sigma2[1] = sigma2_0;
  Mu[1] = mu1;
  Mu[2] = mu2;

  mu1 ~ normal(0.5, 10000);
  mu2 ~ normal(0.5, 10000);

  a1_raw ~ beta(1, 1);
  b1_raw ~ beta(1, 1);
  a2_raw ~ beta(1, 1);
  b2_raw ~ beta(1, 1);

  w1 ~ normal(1.0, 10000);
  w2 ~ normal(1.0, 10000);
  r ~ beta(1, 1);

  for (t in 2:T) {
    sigma1[t] = sqrt(
      w1
      + a1 * square(Y[t - 1][1] - mu1)
      + b1 * square(sigma1[t - 1])
    );

    sigma2[t] = sqrt(
      w2
      + a2 * square(Y[t - 1][2] - mu2)
      + b2 * square(sigma2[t - 1])
    );

    Sigma[1, 1] = square(sigma1[t]);
    Sigma[2, 2] = square(sigma2[t]);
    Sigma[1, 2] = r * sigma1[t] * sigma2[t];
    Sigma[2, 1] = Sigma[1, 2];

    target += power * multi_normal_lpdf(Y[t] | Mu, Sigma);
  }
}
