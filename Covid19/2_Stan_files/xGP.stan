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

	int<lower=0> deaths[T, N_age_groups];
	matrix<lower=0>[N_age_groups, N_age_groups] Cont_mat;
	
	vector<lower=0>[N_age_groups] pop;
	
	vector<lower=0>[T] I_D; // Discretized infection to death distribution.  
	vector<lower=0>[T] GT;  // Discretized generation time interval.
	real<lower=0,upper=1> ifr[N_age_groups];
	
	int<lower=0> epi_start;
	
	real dI; 
}

transformed data {
	real time[T] = linspaced_array(T, 1, T);
	int M = (T + gp_days_cp - 1) / gp_days_cp;  // ceil(T/gp_days_cp), number of blocks
	real time_adjusted[M]= linspaced_array(M, 1, M);

	int time_blocks[T];
	for (t in 1:T) time_blocks[t] = (t - 1) / gp_days_cp + 1;  // integer division groups days

	row_vector<lower = 0>[T] I_D_rev; // Reversed discretized infection-to-death distribution
	row_vector<lower = 0>[T] GT_rev;  // Reversed discretized generation time interval

	for(i in 1:T) {
	  I_D_rev[i] = I_D[T - i + 1];
	  GT_rev[i]  = GT[T - i + 1];
	}// End for
}

parameters {
	// matrix<lower=0>[N_age_groups, N0] cases0;
	row_vector<lower=0>[N_age_groups] cases0; // Expected number of daily infections for each group
	real<lower=0> phi;
	real<lower=0> tau;
	
	real eta_age_group[N_age_groups*M];         // Non-daily frequency
	// matrix[N_age_groups, M] eta_age_group;
	vector[M] eta_country;
	
	real<lower=0> sigma_age_group;
	real<lower=0> lambda_age_group;
	real<lower=0> sigma_country;
	real<lower=0> lambda_country;
}

transformed parameters {
  matrix[M, M] K_country_adjusted     = add_diag(cov_exp_quad(time_adjusted, sigma_country,   lambda_country),   1e-6);
  matrix[M, M] K_age_group_adjusted   = add_diag(cov_exp_quad(time_adjusted, sigma_age_group, lambda_age_group), 1e-6);
  matrix[M, N_age_groups] eta_age_mat = to_matrix(eta_age_group, M, N_age_groups);

  matrix[M, M] L_K_country_adjusted;
  matrix[M, M] L_K_age_group;
  
  matrix[M, N_age_groups] logit_p_adjusted;				  
  matrix[T, N_age_groups] logit_p;
  
  row_vector[N_age_groups] conv_tmp;

  matrix<lower=0, upper=1>[T, N_age_groups] p;
  matrix<lower=0>[T, N_age_groups] cases;
  matrix<lower=0>[T, N_age_groups] E_deaths;
  
  L_K_country_adjusted = cholesky_decompose(K_country_adjusted);
  L_K_age_group        = cholesky_decompose(K_age_group_adjusted);

  for (i in 1:N_age_groups) {
    //logit_p_adjusted[, i] = (L_K_age_group * eta_age_mat[,i]'+ L_K_country_adjusted * eta_country)';
    logit_p_adjusted[, i] = L_K_age_group * eta_age_mat[, i] + L_K_country_adjusted * eta_country;
	for (t in 1:T) logit_p[t, i] = logit_p_adjusted[time_blocks[t], i];
  }
  
  p = inv_logit(logit_p);
  
  //---- Initiate expected infections by age in the first N0 days:
  cases[1:N0, ] = rep_matrix(cases0, N0);
  
  for (t in (N0+1):T) {
    for (j in 1:N_age_groups) {
      for (k in 1:N_age_groups) conv_tmp[k] = Cont_mat[k,j] * dI * dot_product( cases[1:t-1,k], tail(GT_rev, t-1) );
	  
	  cases[t,j] = fmax(1e-03, 
			           (1.0 - (rep_row_vector(1.0, t-1) * cases[1:(t-1),j] / pop[j]) ) * p[t,j] * sum(conv_tmp) );
    }
  }
  
  for (i in 1:N_age_groups) E_deaths[1, i] = ifr[i] * cases[1, i] * I_D[T];
  
  for (t in 2:T) for (j in 1:N_age_groups) E_deaths[t,j] = ifr[j] * dot_product( cases[1:t-1, j], tail(I_D_rev, t-1) );

}

model{
  
  target += std_normal_lpdf(eta_country);
  target += std_normal_lpdf(sigma_age_group);
  target += std_normal_lpdf(lambda_age_group);
  target += std_normal_lpdf(sigma_country);
  target += std_normal_lpdf(lambda_country);
  target += exponential_lpdf(tau | 0.03);
  target += std_normal_lpdf(phi);
  
  target += normal_lpdf(eta_age_group | 0, 1);             // eta_noise ~ std_normal();
  target += exponential_lpdf(cases0 | 1 / tau);

  for (i in 1:N_age_groups) {
	  //target += std_normal_lpdf(eta_age_group[i, 1:M]);
	  //target += exponential_lpdf(cases0[i, 1:N0] | 1 / tau);
      target += neg_binomial_2_lpmf( deaths[epi_start:T, i] | E_deaths[epi_start:T, i], phi);
  }
}

generated quantities {
  
  matrix[T - epi_start + 1, N_age_groups] log_like_age; 
  vector[T - epi_start + 1] log_lik; // Log-likelihood vector for loo package.
  
  real deviance;
  
  for ( i in 1:(T - epi_start + 1) ) {
    for (j in 1:N_age_groups) log_like_age[i,j] = neg_binomial_2_lpmf(deaths[epi_start + i - 1, j] | E_deaths[epi_start + i - 1, j], phi);
    log_lik[i] = sum(log_like_age[i,]);
  }// End for
  
  //---- Deviance:
  deviance = (-2) * sum(log_lik);
 
}
