functions {
  vector matrix_to_vector(matrix X, int num_rows, int num_cols) {
    return to_vector(to_row_vector(X)');
  }
}

data {
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
}

transformed data {
  real time[T] = linspaced_array(T, 1, T);
}

parameters {
  real<lower=0> sigma_age_group;
  real<lower=0> lambda_age_group;
  real<lower=0> sigma_country;
  real<lower=0> lambda_country;
  matrix[N_age_groups, T] eta_age_group;
  vector[T] eta_country;
  matrix<lower=0>[N_age_groups, N0] cases0;
  real<lower=0> phi;
  real<lower=0> tau;
}

transformed parameters {
  matrix[T, T] K_country = add_diag(cov_exp_quad(time, sigma_country, lambda_country), 1e-6);
  matrix[T, T] K_age_group= add_diag(cov_exp_quad(time, sigma_age_group, lambda_age_group), 1e-6);
  matrix[T, T] L_K_age_group;
  matrix[T, T] L_K_country;

  matrix[N_age_groups, T] logit_p;
  matrix<lower=0, upper=1>[N_age_groups, T] p;
  matrix<lower=0>[N_age_groups, T] cases;
  matrix<lower=0>[N_age_groups, T] E_deaths;
  
  L_K_age_group=cholesky_decompose(K_age_group);
  L_K_country=cholesky_decompose(K_country);
  for (i in 1:N_age_groups) {
    logit_p[i] = eta_age_group[i] * L_K_age_group +to_row_vector(eta_country)*L_K_country;
  }
  
  p = inv_logit(logit_p);
  cases[, 1:N0] = cases0;

  for (t in (N0+1):T) {
    for (i in 1:N_age_groups) {
      vector[N_age_groups] conv_tmp = p[, t] .* Cont_mat[, i] .* rows_dot_product(cases[, 1:(t-1)], rep_matrix(to_row_vector(s_int_rev1[(T-t+2):T]),N_age_groups));
      cases[i, t] = (1 - sum(cases[i, 1:(t-1)]) / pop[i]) * sum(conv_tmp);
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

model {
for(i in 1:N_age_groups){
  eta_age_group[i,1:T]~std_normal();
  cases0[i,1:N0] ~ exponential(1/tau);
  //likelihood
  deaths[i,epi_start:T] ~ neg_binomial_2(E_deaths[i,epi_start:T], phi);

}
  eta_country~std_normal();
  sigma_age_group ~std_normal();
  lambda_age_group ~ std_normal();
  sigma_country ~ std_normal();
  lambda_country ~ std_normal();
  tau ~ exponential(0.03);
  phi ~ std_normal();
}

generated quantities {
  real Deviance;
  vector[N_age_groups * (T - epi_start + 1)] log_likelihood;
  matrix[N_age_groups, (T - epi_start + 1)] log_likelihood_matrix;
  matrix<lower=0>[N_age_groups, T] Rt;
  real total_E_deaths[T];
  real total_cases[T];

  
  for (t in 1:(T - epi_start + 1)) {
    for (i in 1:N_age_groups) {
      log_likelihood_matrix[i, t] = neg_binomial_2_lpmf(deaths[i, epi_start + t - 1] | E_deaths[i, epi_start + t - 1], phi);
    }
  }
  log_likelihood = matrix_to_vector(log_likelihood_matrix, N_age_groups, (T - epi_start + 1));

  for (t in 1:T) {
    for (i in 1:N_age_groups) {
      vector[N_age_groups] Rt_tmp = (1 - rows_dot_product(cases[, 1:(t-1)], rep_matrix(1.0,N_age_groups,t-1))./ pop) .* p[, t] .* to_vector(Cont_mat[i]);
      Rt[i, t] = sum(Rt_tmp);
    }
    total_E_deaths[t]=sum(E_deaths[1:N_age_groups,t]);
    total_cases[t]=sum(cases[1:N_age_groups,t]);

  }

  Deviance = -2 * sum(log_likelihood);
}
