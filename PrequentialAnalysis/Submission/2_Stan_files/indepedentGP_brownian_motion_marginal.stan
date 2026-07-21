data{
  int<lower=0> N; //number of timepoints-obervations 

  vector[N] Y; // the observed directory 

  real Times[N]; //Different Timepoints
  
  
}
transformed data {

  vector[N] mu = rep_vector(0, N);


  real delta = 1e-9;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> sigma_e;
}
transformed parameters{
  matrix[N, N] L_K;
  matrix[N, N] K;
  // Build Brownian motion covariance matrix: K[i,j] = sigma^2 * min(t_i, t_j)
  for (i in 1:N) {
    for (j in 1:N) {
      K[i, j] = square(sigma) * fmin(Times[i], Times[j]);
    }
  }
  L_K= cholesky_decompose(K+diag_matrix(rep_vector(square(sigma_e) + delta, N)));
  
}
model {
  //priors
  sigma~std_normal();
  sigma_e~std_normal();
  //marginal likelihood
  Y ~ multi_normal_cholesky(mu, L_K);
}

