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
	
	/*
	matrix[W,C] weather[K];        // Weather covariates - Array of K matrices (islands), each having W rows and C columns

	//array[K] matrix[W,C] weather;  // Weather covariates - Array of K matrices (islands), each having W rows and C columns
	// For array indexing between R and Stan, see https://discourse.mc-stan.org/t/dimension-mismatch-when-passing-array-of-matrix-from-r-to-stan/22404/7
    */
	
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
	
	//---- Debugging:
	int doprint;
}

parameters {
	// Logit of the distribution of rho for CHIKV
	//real logitrhoc_raw[K];   // Island-specific logit reporting rate for CHIKV
	//real logitrhoc[K];   // Island-specific logit reporting rate for CHIKV

	//real mu_rhoc;             
	//real<lower = 0> sigma_rhoc;
	
	// Random parameters
    real mu_b0c;     
    real<lower=0> sigma_b0c;
    real b0ic[K]; // island-specific base transmission for CHIKV

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
    matrix[C, K] bW_mat = to_matrix(bW, C, K);

	for(k in 1:K) { // Across islands
	
		////logitrhoc[k] = mu_rhoc + sigma_rhoc * logitrhoc_raw[k];     // Non-centered parameterisation
		//rhoc[k]      = exp(logitrhoc[k]) / (1 + exp(logitrhoc[k]));

	    theta[,k]    = weather[,(C*k - C + 1):(C*k)] * bW_mat[,k];
				
		for (i in 1:W){
			beta[i,k] = exp( b0ic[k] + theta[i,k] );
			lp[i,k]   = beta[i,k] * Ostar_t[i,k] * ( 1 - sumO_t[i,k] / (pop[k]) ) + 0.001; //(rhoc[k] * pop[k])
		}// End for
	}// End for
	
	//---- Debug:
	if(doprint != 0) {
	  print("phi: ", phi);
	 // print("rhoc: ", rhoc);
	  print("theta: ", theta); 
	  print("beta: ", beta); 
	  print("lp: ", lp); 
	}// End debugging
}

model {
	
	//---- Hyperparameter priors 
	//mu_rhoc    ~ std_normal();
	//sigma_rhoc ~ inv_gamma(10,10);
	
	target  += student_t_lpdf(mu_b0c | 5, 0, 2.5);
	target  += cauchy_lpdf(sigma_b0c | 0, 2.5);
    target  += normal_lpdf(b0ic | mu_b0c, sigma_b0c);

	//---- (Rest of) Random parameter priors
	//logitrhoc_raw ~ std_normal();
	//logitrhoc ~ normal(mu_rhoc, sigma_rhoc);
	
	//---- Fixed parameters priors
	target  += student_t_lpdf(bW | 5,0,2.5);
	
	//---- Overdispersion prior
	if (prior_dist_nb_dispersion == 1)      target  += normal_lpdf(phi | prior_mean_nb_dispersion, prior_scale_nb_dispersion);            // Half-Normal
	else if (prior_dist_nb_dispersion == 2) target  += cauchy_lpdf(phi | prior_mean_nb_dispersion, prior_scale_nb_dispersion);                            // Half-Cauchy
	else if (prior_dist_nb_dispersion == 3) target  += student_t_lpdf(phi | prior_df_nb_dispersion, prior_mean_nb_dispersion, prior_scale_nb_dispersion); // Half Student-t
	else if (prior_dist_nb_dispersion == 4) target  += gamma_lpdf(phi | prior_shape_nb_dispersion, prior_rate_nb_dispersion);                             // Gamma
	else if (prior_dist_nb_dispersion == 5) target  += exponential_lpdf(phi | prior_rate_nb_dispersion);                                                  // Exponential

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
	
	for (t in 1:W) {
		for (k in 1:K) {
			if (likelihood_variance_type == 0)      log_like_island[t,k] = neg_binomial_2_lpmf(O_t[t,k] | lp[t,k], phi);
			else if (likelihood_variance_type == 1) log_like_island[t,k] = neg_binomial_2_lpmf(O_t[t,k] | lp[t,k], lp[t,k]/phi);
		}// End for

		log_lik[t] = sum(log_like_island[t,]);
	}// End for

	//---- Deviance:
	deviance = (-2) * sum(log_lik);
}
