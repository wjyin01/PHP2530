// dm_regression.stan
data {
  int<lower=1>  N;               // number of samples
  int<lower=1>  K;               // number of taxa/categories
  int<lower=1>  P;               // number of predictors (incl. intercept)
  array[N, K] int<lower=0> Y;    // count matrix (rows = samples)
  matrix[N, P]  X;               // design matrix
}

parameters {
  matrix[P, K]  beta;            // regression coefficients
  vector<lower=0>[N] phi;        // concentration parameter (per sample)
}
transformed parameters {
  matrix[N, K] alpha;            // Dirichlet parameters for each sample
  for (i in 1:N) {
    row_vector[K] eta_row = X[i] * beta;         // row_vector result
    vector[K]     eta     = eta_row';            // convert to column vector
    vector[K]     theta   = softmax(eta);        // simplex
    alpha[i]              = (theta * phi[i])';   // scale and transpose
  }
}
model {
  // Priors
  to_vector(beta) ~ normal(0, 2);       // weakly informative prior
  phi             ~ gamma(2, 0.1);      // prior over concentration

  // Likelihood: Dirichlet-Multinomial
  for (i in 1:N)
    Y[i] ~ dirichlet_multinomial(alpha[i]');
}
generated quantities {
  // Posterior mean compositions
  matrix[N, K] theta_hat;
  for (i in 1:N) {
    row_vector[K] eta_row = X[i] * beta;
    vector[K]     eta     = eta_row';
    vector[K]     theta   = softmax(eta);
    theta_hat[i]          = theta';   // transpose to row_vector
  }
}
