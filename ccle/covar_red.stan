// generated with brms 2.17.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += lognormal_lpdf(b | 0, 1);
  lprior += gamma_lpdf(shape | 0.01, 0.01);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = X * b;
    target += neg_binomial_2_lpmf(Y | mu, shape);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
}
