functions {
  vector matrix_to_vector(matrix X, int num_rows, int num_cols) {
    return to_vector(to_row_vector(X)');
  }
}

data {
  int<lower=0> gp_days_cp;
  int<lower=0> N0;
  int<lower=1> T;
  int<lower=1> N_age_groups;
  int<lower=0> deaths[N_age_groups, T];
  matrix<lower=0>[N_age_groups, N_age_groups] Cont_mat;
  vector<lower=0>[N_age_groups] pop;
  vector<lower=0>[T] s_int_rev1;
  vector<lower=0>[T] s_int_rev2;
  real<lower=0,upper=1> ifr[N_age_groups];
  int<lower=0> epi_start;
  real dI; 
}

transformed data {
  real time[T] = linspaced_array(T, 1, T);
  int M = (T + gp_days_cp - 1) / gp_days_cp;  // ceil(T/gp_days_cp), number of blocks
  real time_adjusted[M]= linspaced_array(M, 1, M);
  
  int time_blocks[T];
  for (t in 1:T) {
    time_blocks[t] = (t - 1) / gp_days_cp + 1;  // integer division groups days
  }
}

parameters {
  vector<lower=0>[N_age_groups] sigma_age_group;
  vector<lower=0>[N_age_groups] lambda_age_group;
  real<lower=0> sigma_country;
  real<lower=0> lambda_country;
  matrix[N_age_groups, M] eta_age_group;
  vector[M] eta_country;

  matrix<lower=0>[N_age_groups, N0] cases0;
  real<lower=0> phi;
  real<lower=0> tau;
}

transformed parameters {
  matrix[M, M] K_country_adjusted = add_diag(cov_exp_quad(time_adjusted, sigma_country, lambda_country), 1e-6);
  matrix[M, M] K_age_group_adjusted[N_age_groups];
  matrix[N_age_groups, M] logit_p_adjusted;
  matrix[N_age_groups, T] logit_p;
  matrix<lower=0, upper=1>[N_age_groups, T] p;
  matrix<lower=0>[N_age_groups, T] cases;
  matrix<lower=0>[N_age_groups, T] E_deaths;
  matrix[M,M] L_K_country_adjusted;
  
  
  L_K_country_adjusted=cholesky_decompose(K_country_adjusted);

  for (i in 1:N_age_groups) {
    K_age_group_adjusted [i] = add_diag(cov_exp_quad(time_adjusted, sigma_age_group[i], lambda_age_group[i]), 1e-6);
    logit_p_adjusted[i]      = (cholesky_decompose(K_age_group_adjusted[i]) * eta_age_group[i]'+ L_K_country_adjusted * eta_country)';
    for (t in 1:T) logit_p[i,t] = logit_p_adjusted[i,time_blocks[t]];
  }
  
  p = inv_logit(logit_p);
  cases[, 1:N0] = cases0;

  for (t in (N0+1):T) {
    for (i in 1:N_age_groups) {
      vector[N_age_groups] conv_tmp = p[, t] .* Cont_mat[, i] .* rows_dot_product(cases[, 1:(t-1)], rep_matrix(to_row_vector(s_int_rev1[(T-t+2):T]),N_age_groups));
      cases[i, t] = fmax(1e-03,
	                    (1 - sum(cases[i, 1:(t-1)]) / pop[i]) * dI * sum(conv_tmp)
						);
    }
  }
  
  for (i in 1:N_age_groups) {
    E_deaths[i, 1] = ifr[i] * cases[i, 1] * s_int_rev2[T];
  }
  for (t in 2:T) {
    for (i in 1:N_age_groups) {
      E_deaths[i, t] = ifr[i] * (cases[i, 1:(t-1)] * s_int_rev2[(T-t+2):T]);
    }
  }
}


model{
  for (i in 1:N_age_groups) {
    target += std_normal_lpdf(eta_age_group[i, 1:M]);
    target += exponential_lpdf(cases0[i, 1:N0] | 1 / tau);
    //likelihood
    target += neg_binomial_2_lpmf(deaths[i, epi_start:T] | E_deaths[i, epi_start:T], phi);
  }

  target += std_normal_lpdf(eta_country);
  target += std_normal_lpdf(sigma_age_group);
  target += std_normal_lpdf(lambda_age_group);
  target += std_normal_lpdf(sigma_country);
  target += std_normal_lpdf(lambda_country);
  target += exponential_lpdf(tau | 0.03);
  target += std_normal_lpdf(phi);
}

generated quantities {
  matrix[T - epi_start + 1, N_age_groups] log_like_age; 
  vector[T - epi_start + 1] log_lik; // Log-likelihood vector for loo package.
  
  real deviance;
  
  for ( i in 1:(T - epi_start + 1) ) {
    for (j in 1:N_age_groups) log_like_age[i,j] = neg_binomial_2_lpmf(deaths[j, epi_start + i - 1] | E_deaths[j, epi_start + i - 1], phi);
    log_lik[i] = sum(log_like_age[i,]);
  }// End for
  
  //---- Deviance:
  deviance = (-2) * sum(log_lik);
}
