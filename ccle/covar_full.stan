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
  lprior += normal_lpdf(b[1] | 0, 0.5);
  lprior += lognormal_lpdf(b[2] | 0, 1);
  lprior += lognormal_lpdf(b[3] | 0, 1);
  lprior += lognormal_lpdf(b[4] | 0, 1);
  lprior += lognormal_lpdf(b[5] | 0, 1);
  lprior += lognormal_lpdf(b[6] | 0, 1);
  lprior += lognormal_lpdf(b[7] | 0, 1);
  lprior += lognormal_lpdf(b[8] | 0, 1);
  lprior += lognormal_lpdf(b[9] | 0, 1);
  lprior += lognormal_lpdf(b[10] | 0, 1);
  lprior += lognormal_lpdf(b[11] | 0, 1);
  lprior += lognormal_lpdf(b[12] | 0, 1);
  lprior += lognormal_lpdf(b[13] | 0, 1);
  lprior += lognormal_lpdf(b[14] | 0, 1);
  lprior += lognormal_lpdf(b[15] | 0, 1);
  lprior += lognormal_lpdf(b[16] | 0, 1);
  lprior += lognormal_lpdf(b[17] | 0, 1);
  lprior += lognormal_lpdf(b[18] | 0, 1);
  lprior += lognormal_lpdf(b[19] | 0, 1);
  lprior += lognormal_lpdf(b[20] | 0, 1);
  lprior += lognormal_lpdf(b[21] | 0, 1);
  lprior += lognormal_lpdf(b[22] | 0, 1);
  lprior += lognormal_lpdf(b[23] | 0, 1);
  lprior += lognormal_lpdf(b[24] | 0, 1);
  lprior += lognormal_lpdf(b[25] | 0, 1);
  lprior += lognormal_lpdf(b[26] | 0, 1);
  lprior += lognormal_lpdf(b[27] | 0, 1);
  lprior += lognormal_lpdf(b[28] | 0, 1);
  lprior += lognormal_lpdf(b[29] | 0, 1);
  lprior += lognormal_lpdf(b[30] | 0, 1);
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
