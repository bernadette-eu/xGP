functions {
  // Removed user-defined SE kernel function
}

data {
  int inference;
  int<lower=1> W; //number of records
  int<lower=1> K; //number of islands
  int<lower=0> O_t[W, K];
  matrix<lower=0>[W, K] Ostar_t;
  matrix<lower=0>[W, K] sumO_t;
  int<lower=0> pop[K];
  int<lower=0> C;
  matrix[W,C*K] weather;
  real<lower=0> timepoints[W];
  
  	/*
	//---- Negative binomial variance formulation:
	0 = variance as a quadratic function of the mean
	1 = variance as a linear function of the mean
	*/
	int<lower = 0, upper = 1> likelihood_variance_type;

	/*
	//---- Prior hyperparameters for NegativeBinomial dispersion parameter:
	1 = Half-Normal(prior_mean_nb_dispersion, prior_scale_nb_dispersion)
	2 = Half-Cauchy(prior_mean_nb_dispersion, prior_scale_nb_dispersion)
	3 = Half Student-t(prior_df_nb_dispersion, prior_mean_nb_dispersion, prior_scale_nb_dispersion)
	4 = Gamma(prior_shape_nb_dispersion, prior_rate_nb_dispersion)
	5 = Exponential(prior_rate_nb_dispersion)
	*/
	int<lower = 1, upper = 5> prior_dist_nb_dispersion;
	real<lower=0> prior_mean_nb_dispersion;
	real<lower=0> prior_scale_nb_dispersion;
	real<lower=0> prior_df_nb_dispersion;
	real<lower=0> prior_shape_nb_dispersion;
	real<lower=0> prior_rate_nb_dispersion;

  	/*
	//---- Prior hyperparameters for volatility parameters:
	1 = Half-Normal(prior_mean_volatility, prior_scale_volatility)
	2 = Half-Cauchy(prior_mean_volatility, prior_scale_volatility)
	3 = Half Student-t(prior_df_volatility, prior_mean_volatility, prior_scale_volatility)
	4 = Gamma(prior_shape_volatility, prior_rate_volatility)
	5 = Exponential(prior_rate_volatility)
	*/
	int<lower = 1, upper = 5> prior_dist_sigma;
	real<lower=0> prior_mean_sigma;
	real<lower=0> prior_scale_sigma;
	real<lower=0> prior_df_sigma;
	real<lower=0> prior_shape_sigma;
	real<lower=0> prior_rate_sigma;
}

parameters {
  matrix[W, K] x_gp_raw;
  real<lower=0> sigma[K];
  real<lower=0> length_scale[K];
  real bW[C*K];
  real<lower = 0> phi;
}

transformed parameters {
  matrix[W, K] theta;
  matrix<lower = 0>[W, K] beta;
  matrix<lower = 0>[W, K] lp;
  matrix[C, K] bW_mat = to_matrix(bW, C, K);
  matrix[W, K] x_trajectory;

  for(k in 1:K) {
    matrix[W, W] K_se = add_diag(cov_exp_quad(timepoints, sigma[k], length_scale[k]), 1e-6);
    matrix[W, W] L_se = cholesky_decompose(K_se);
    x_trajectory[,k] = L_se * x_gp_raw[,k];
  }

  for(k in 1:K) {
    theta[,k] = weather[,(C*k - C + 1):C*k] * bW_mat[,k];
    for (i in 1:W){
      beta[i,k] = exp(x_trajectory[i,k] + theta[i,k]);
      lp[i,k] = beta[i,k] * Ostar_t[i,k] * (1 - sumO_t[i,k] / (pop[k])) + 0.001;
    }
  }
}

model {
  
	target  += normal_lpdf(to_vector(x_gp_raw) | 0, 1);
  
	if (prior_dist_sigma == 1)       target  += normal_lpdf(sigma | prior_mean_sigma, prior_scale_sigma);                    // Half-Normal
	else if (prior_dist_sigma == 2)  target  += cauchy_lpdf(sigma | prior_mean_sigma, prior_scale_sigma);                    // Half-Cauchy
	else if (prior_dist_sigma == 3)  target  += student_t_lpdf(sigma | prior_df_sigma, prior_mean_sigma, prior_scale_sigma); // Half Student-t
	else if (prior_dist_sigma == 4)  target  += gamma_lpdf(sigma | prior_shape_sigma, prior_rate_sigma);                     // Gamma
	else if (prior_dist_sigma == 5)  target  += exponential_lpdf(sigma | prior_rate_sigma);                                  // Exponential

	for (k in 1:K) target  += normal_lpdf(length_scale[k] | 0, 1); //length_scale[k] ~ std_normal();
	target  += student_t_lpdf(bW | 5, 0, 2.5);
  
	if (prior_dist_nb_dispersion == 1)      target  += normal_lpdf(phi | prior_mean_nb_dispersion, prior_scale_nb_dispersion);            // Half-Normal
	else if (prior_dist_nb_dispersion == 2) target  += cauchy_lpdf(phi | prior_mean_nb_dispersion, prior_scale_nb_dispersion);                            // Half-Cauchy
	else if (prior_dist_nb_dispersion == 3) target  += student_t_lpdf(phi | prior_df_nb_dispersion, prior_mean_nb_dispersion, prior_scale_nb_dispersion); // Half Student-t
	else if (prior_dist_nb_dispersion == 4) target  += gamma_lpdf(phi | prior_shape_nb_dispersion, prior_rate_nb_dispersion);                             // Gamma
	else if (prior_dist_nb_dispersion == 5) target  += exponential_lpdf(phi | prior_rate_nb_dispersion);                                                  // Exponential

	if(inference == 1){
		for(t in 1:W) {
		  for (k in 1:K) {
			if (likelihood_variance_type == 0)      target += neg_binomial_2_lpmf(O_t[t,k] | lp[t,k], phi);
			else if (likelihood_variance_type == 1) target += neg_binomial_2_lpmf(O_t[t,k] | lp[t,k], lp[t,k]/phi);
		  }
		}
  }
}

generated quantities {
  real deviance;
  matrix[W,K] log_like_island;
  vector[W] log_lik;
  matrix[W,K] rt_eff_island;
  matrix[W,K] rt_eff_island2;
	
	for (k in 1:K) rt_eff_island2[1,k] = (1-lp[1,k]/pop[k])*exp(x_trajectory[1,k]);

	for (t in 1:W) {
		for (k in 1:K) {
			if (likelihood_variance_type == 0)      log_like_island[t,k] = neg_binomial_2_lpmf(O_t[t,k] | lp[t,k], phi);
			else if (likelihood_variance_type == 1) log_like_island[t,k] = neg_binomial_2_lpmf(O_t[t,k] | lp[t,k], lp[t,k]/phi);
		    
			rt_eff_island[t,k]  = (1-lp[t,k]/pop[k])*exp(x_trajectory[t,k]);
			if (t >=2) rt_eff_island2[t,k] = (1-rep_row_vector(1.0, t-1) * lp[1:(t-1),k]/pop[k])*exp(x_trajectory[t,k]);

		}// End for

		log_lik[t] = sum(log_like_island[t,]);
	}// End for

	//---- Deviance:
	deviance = (-2) * sum(log_lik);
}
