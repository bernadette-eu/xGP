functions {
		
	//---- Formulation of the covariance matrix for \Delta A_t:
	matrix cov_mat_est(int A,
					   vector country_vars
					   ){

		matrix[A, A] y; 
		matrix[A, A] country_diag = diag_matrix(country_vars);

		y = country_diag;

		return y;
			
	}// End function

}

data {
	int inference;
	int<lower=1> W;                // Number of records
	int<lower=1> K;                // Number of islands
	//int<lower = 1, upper = 7> ecr_changes;
	//int n_changes;
	//int n_remainder;
	
	//int<lower=1> island[W];     // Identification of island 1=TAHITI 2=TUAMOTU 3=MOOREA 4=SLV 5=MARQUISES 6=AUSTRALES 7=SAINT MARTIN 8=MARTINIQUE 9=GUADELOUPE

	int<lower=0> O_t[W, K];        // Number of reported cases at time t
	matrix<lower=0>[W, K] Ostar_t; // Exposure at time t
	matrix<lower=0>[W, K] sumO_t;  // Cumulative number of reported cases at time t

	int<lower=0> pop[K];           // Island population
	int<lower=0> C;                // Number of weather covariables
	
	matrix[W,C*K] weather;         // Weather covariates - Array of K matrices (islands), each having W rows and C columns

	//matrix[W,C] weather[K];        // Weather covariates - Array of K matrices (islands), each having W rows and C columns
	// For array indexing between R and Stan, see https://discourse.mc-stan.org/t/dimension-mismatch-when-passing-array-of-matrix-from-r-to-stan/22404/7

	/*
	//---- Negative binomial variance formulation:
	0 = variance as a quadratic function of the mean
	1 = variance as a linear function of the mean
	*/
	int<lower = 0, upper = 1> likelihood_variance_type;

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
}

parameters {
	// Logit of the distribution of rho for CHIKV
	//real logitrhoc_raw[K];   // Island-specific logit reporting rate for CHIKV
	//real mu_rhoc;             
	//real<lower = 0> sigma_rhoc;
	
	// Random parameters
	real x_noise[W*K];
	real<lower = 0> sigma[K];         // Volatility of each island

	// Fixed parameters (weather)
	real bW[C*K]; // Effects of the C weather covariates for each island

	// Dispersion parameter
	real<lower = 0> phi;
}

transformed parameters {
	//real logitrhoc[K];                 // Island-specific logit reporting rate for CHIKV
	//real<lower = 0,upper = 1> rhoc[K]; // reporting rate by island
	matrix[W, K] theta;
	
	matrix<lower = 0>[W, K] beta; // Transmission
	matrix<lower = 0>[W, K] lp;   // Expexcted incidence in time per island

    vector[K] island_vars;
	matrix[W, K] delta_x_mat;  
	matrix[K, K] cov_mat; 
	
	matrix[W, K] x_noise_mat  = to_matrix(x_noise, W, K);
	matrix[W, K] x_trajectory = rep_matrix(0.0, W, K); // Island-specific base transmission for CHIKV
	matrix[C, K] bW_mat       = to_matrix(bW, C, K);

	//---- Cholesky decomposition of the covariance matrix.
	for (j in 1:K) island_vars[j] = sigma[j]^2;
	cov_mat = cov_mat_est(K, island_vars);
	
	// Non-centered parameterisation:
	// https://occasionaldivergences.com/posts/non-centered/
	for(t in 1:W) delta_x_mat[t,] = to_row_vector( cholesky_decompose(cov_mat) * x_noise_mat[t,]'); // x_noise_mat[t,] * cov_mat; 

	//---- Affine transformations for xGBMs (Non-central parameterisation):
	// Initiate x_mat at t = 1, assuming that x_{0,j} = 0, j = 1, \ldots, K:

	for(k in 1:K) { // Across islands
	
		//logitrhoc[k] = mu_rhoc + sigma_rhoc * logitrhoc_raw[k];     // Non-centered parameterisation
		//rhoc[k]      = exp(logitrhoc[k]) / (1 + exp(logitrhoc[k]));

	    theta[,k]    = weather[,(C*k - C + 1):C*k] * bW_mat[,k];

		// Island-specific base transmission  
		x_trajectory[1,k]   = delta_x_mat[1,k];
		x_trajectory[2:W,k] = x_trajectory[1:(W - 1),k] + delta_x_mat[2:W,k];
				
		for (i in 1:W){
			beta[i,k] = exp( x_trajectory[i,k] + theta[i,k] );
			lp[i,k]   = beta[i,k] * Ostar_t[i,k] * ( 1 - sumO_t[i,k] / (pop[k]) ) + 0.001; //(rhoc[k] * pop[k])
		}// End for
	}// End for
}

model {
	
	//---- Hyperparameter priors 
	//mu_rhoc    ~ std_normal();
	//sigma_rhoc ~ inv_gamma(10,10);
	
	//---- GBM priors:
	x_noise ~ std_normal();

	if (prior_dist_sigma == 1)       sigma ~ normal(prior_mean_sigma, prior_scale_sigma);                    // Half-Normal
	else if (prior_dist_sigma == 2)  sigma ~ cauchy(prior_mean_sigma, prior_scale_sigma);                    // Half-Cauchy
	else if (prior_dist_sigma == 3)  sigma ~ student_t(prior_df_sigma, prior_mean_sigma, prior_scale_sigma); // Half Student-t
	else if (prior_dist_sigma == 4)  sigma ~ gamma(prior_shape_sigma, prior_rate_sigma);                     // Gamma
	else if (prior_dist_sigma == 5)  sigma ~ exponential(prior_rate_sigma);                                  // Exponential                                 // Exponential

	//---- (Rest of) Random parameter priors
	//logitrhoc_raw ~ std_normal();
	
	//---- Fixed parameters priors
	bW ~ student_t(5,0,2.5);
	
	//---- Overdispersion prior
	if (prior_dist_nb_dispersion == 1)      phi ~ normal(prior_mean_nb_dispersion, prior_scale_nb_dispersion);                            // Half-Normal
	else if (prior_dist_nb_dispersion == 2) phi ~ cauchy(prior_mean_nb_dispersion, prior_scale_nb_dispersion);                            // Half-Cauchy
	else if (prior_dist_nb_dispersion == 3) phi ~ student_t(prior_df_nb_dispersion, prior_mean_nb_dispersion, prior_scale_nb_dispersion); // Half Student-t
	else if (prior_dist_nb_dispersion == 4) phi ~ gamma(prior_shape_nb_dispersion, prior_rate_nb_dispersion);                             // Gamma
	else if (prior_dist_nb_dispersion == 5) phi ~ exponential(prior_rate_nb_dispersion);                                                  // Exponential

	if(inference == 1){
		//---- Likelihood:
		for(t in 1:W) {
			for (k in 1:K) {
				if (likelihood_variance_type == 0)      target += neg_binomial_2_lpmf( O_t[t,k] | lp[t,k], phi);
				else if (likelihood_variance_type == 1) target += neg_binomial_2_lpmf( O_t[t,k] | lp[t,k], lp[t,k]/phi);
			}// End for
		}// End for
	}// End if
	
}

generated quantities {
	real deviance;     // Deviance
	matrix[W,K] log_like_island;
	vector[W] log_lik; // Log-likelihood vector for LOOIC computation
	matrix[W,K] rt_eff_island;
	matrix[W,K] rt_eff_island2;
	
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
