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
  real<lower=0> length_scale;
  real<lower=0> sigma;
  real<lower=0> sigma_e;
}
transformed parameters{
  matrix[N, N] L_K;
  matrix[N, N] K;
  K=gp_exp_quad_cov(Times, sigma, length_scale);
  L_K= cholesky_decompose(K+diag_matrix(rep_vector(square(sigma_e) + delta, N)));
  
}
model {
  //priors
  length_scale~std_normal();
  sigma~std_normal();
  sigma_e~std_normal();
  //marginal likelihood
  Y ~ multi_normal_cholesky(mu, L_K);
}

