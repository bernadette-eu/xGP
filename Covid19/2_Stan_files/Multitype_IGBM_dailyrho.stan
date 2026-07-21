functions {
	vector matrix_to_vector(matrix X, int num_rows, int num_cols) {
		return to_vector(to_row_vector(X)');
	}
  
	// Source: https://stackoverflow.com/questions/60810640/what-is-the-equivalent-to-the-r-function-repx-each-n-in-stan
	real[] rep_each_vector(real[] x, int K) {
		int N = size(x);
		real y[N  *  K];
		int pos = 1;
		
		for (n in 1:N) {
		  for (k in 1:K) {
			y[pos] = x[n];
			pos += 1;
		  }
		}
		return y;
	}// End function
	
	// Source: https://stackoverflow.com/questions/60810640/what-is-the-equivalent-to-the-r-function-repx-each-n-in-stan
	matrix rep_each_matrix(matrix x, int K) {
		int N = rows(x);
	    int M = cols(x);

		matrix[N * K, M] y; //Output: extended matrix
		
		int pos = 1;
		
		for (n in 1:N) {
		  for (k in 1:K) {
			  for (m in 1:M) {
			    y[pos,m] = x[n,m];
		      }
		      pos += 1;
		  }
		}
		return y;
	}// End function
	
    //---- Formulation of the covariance matrix for \Delta A_t:
	matrix cov_mat_est(int A,
					  vector country_vars
					  ){
	
	matrix[A, A] y; 
    matrix[A, A] country_diag = diag_matrix(country_vars);

    y = country_diag;
		 
    return y;
					
	}// End function
				 
  vector rowSum(matrix M){
    int nrow; 
    vector[rows(M)] sums; 
    
    nrow = rows(M); 
	
    for(i in 1:nrow) sums[i] = sum(M[i,]); //sums[i] = sum(row(M,i)); 

    return(sums);
  }
}

data {							  
	//---- Structure:								  
	int<lower = 1> A;                   // Number of groups
	int<lower = 1> N;                   // Number of time points of the analysis period
    int<lower = 1> n_weeks;             // Number of weeks
    int gbm_no_days;                    // Possible values: 7, 14, etc. Weekly GBM, biweekly, etc
	
	//---- Data to fit:
	int y_data[N, A];                   // Deaths time series per country
		
	// Structure:								  
	row_vector<lower = 1>[A] n_pop;     // Population per group
	int<lower = 1> N0;                  // Number of days for which to impute infections, N0 < N
	int EpidemicStart;

	vector<lower = 0,upper = 1>[A] ifr; // Infection-fatality rate per country
	vector<lower = 0>[N] I_D;           // Discretized infection to death distribution.     
	vector<lower = 0>[N] GT;            // Discretized generation time interval.          
    matrix[A, A] contact_matrix;

	// Priors:
	real p_phiD;
	real p_sigma_x;

	//---- Debugging:
	int inference; // 0: simulating from priors; 1: fit to data
	int doprint;	
}

transformed data {
	row_vector<lower = 0>[N] I_D_rev; // Reversed discretized infection-to-death distribution
	row_vector<lower = 0>[N] GT_rev;  // Reversed discretized generation time interval
  
	vector[A] ones_vector_A = rep_vector(1.0, A);
	vector[n_weeks] A_vector = rep_vector(A*1.0, n_weeks);
	
	for(i in 1:N) {
	  I_D_rev[i] = I_D[N - i + 1];
	  GT_rev[i]  = GT[N - i + 1];
	}// End for
}

parameters {
    row_vector<lower=0>[A] E_cases_N0; // Expected number of daily infections for each group
	real<lower = 0> phiD;              // Likelihood overdispersion parameter
	real eta_noise[N*A];               // Daily frequency
	real<lower = 0> sigma_x[A];        // Volatility of each country
}

transformed parameters{  
	
	//---- Variables in output:
	matrix[A, A]       cov_mat;    
	matrix[N, A] delta_x_mat;  
	vector[A]          country_vars;

    matrix<lower = 0, upper = 1>[N, A] rho_daily;
	
	matrix<lower = 0>[N, A] E_deathsAge;
	row_vector[A] convolution_tmp;
	
    //---- Transformed parameters for the GBM (Non-central parameterisation):                           			
	matrix[N, A] eta_noise_mat = to_matrix(eta_noise, N, A);
	matrix[N, A] x_mat         = rep_matrix(0.0, N, A);
    matrix[N, A] E_casesAge    = rep_matrix(0.0, N, A );

    //---- Cholesky decomposition of the covariance matrix.
    for (j in 1:A) country_vars[j] = sigma_x[j]^2;
	cov_mat = cov_mat_est(A, country_vars);

	for(t in 1:N) delta_x_mat[t,] = to_row_vector( cholesky_decompose(cov_mat) * eta_noise_mat[t,]');

	//---- Affine transformations for xGBMs:
	// Initiate x_mat at t = 1, assuming that x_{0,j} = 0, j = 1, \ldots, A:
    x_mat[1,:] = delta_x_mat[1,:];
	
	/*
    // This is equivalent to:
	for (j in 1:A) for(t in 2:n_weeks) x_mat[t,j] = delta_x_mat[t,j] + x_mat[t-1,j];
	*/
	for (j in 1:A) x_mat[2:N,j] = delta_x_mat[2:N,j] + x_mat[1:(N-1),j];

    //---- Calculate the probability that a contact with an infectious person leads to infection in group j,
	// by taking the inverse logit of x_mat. Then convert to daily granularity:
	rho_daily = inv_logit(x_mat);
	
	//for (j in 1:A) for(t in 1:N) rho_daily[t,j] = max([rho_daily[t,j], 0]);

	//---- Initiate expected infections by age in the first N0 days:
	E_casesAge[1:N0,] = rep_matrix( E_cases_N0, N0);
    E_deathsAge[1,:]  = rep_row_vector(1e-15, A) .* E_casesAge[1,:];

	//---- Estimation of expected infections:
	for( t in (N0+1):N ){
		//print("t: ", t);

		for (j in 1:A) {
			/*
			print("j: ", j);
		    print("prop_susceptibleByAge, j: ", prop_susceptibleByAge);

			prop_susceptibleByAge = 1.0 - (rep_row_vector(1.0, t-1) * E_casesAge[1:(t-1),j] / n_pop[j]);
            if(prop_susceptibleByAge < 0) prop_susceptibleByAge = 0;
			*/
			
			for (k in 1:A) convolution_tmp[k] = (1.0 - (rep_row_vector(1.0, t-1) * E_casesAge[1:(t-1),j] / n_pop[j]) ) * 
												rho_daily[t,j] * 
												contact_matrix[j,k] * 
												dot_product( E_casesAge[1:t-1,k], tail(GT_rev, t-1) );
			//for (k in 1:A) convolution_tmp[k]  = contact_matrix[j,k] * dot_product( E_casesAge[1:t-1,k], tail(GT_rev, t-1) );
			
			E_casesAge[t,j] = sum(convolution_tmp);
			//E_casesAge[t,j] = prop_susceptibleByAge * rho_daily[t,j] * sum(convolution_tmp);
		}// End for
		
		//print("convolution_tmp: ",  convolution_tmp);
		//print("rho_daily[t,j]: ",  rho_daily[t,]);
	    //print("E_casesAge[t,]: ",  E_casesAge[t,]);

	}// End for

	//---- Estimation of expected deaths:
	for (t in 2:N) for (j in 1:A) E_deathsAge[t,j] = ifr[j] * dot_product( E_casesAge[1:t-1, j], tail(I_D_rev, t-1) );

	//E_deathsAge .*= rep_matrix(ifr, N);
	E_deathsAge += 1e-15;

	// Debug:
	if(doprint != 0) {
	  print("E_casesAge[1:N0,]: ", E_casesAge[1:N0,]);
	  print("rho_daily: ", rho_daily);
	  print("sigma_x: ",   sigma_x);
	  print("Expected Cases [(N0+1):N,]: ",  E_casesAge[(N0+1):N,]);
	  print("phiD: ", phiD);
	  print("Expected Deaths: ", E_deathsAge);
	}// End debugging
}

model {
	
  // Priors:
  E_cases_N0 ~ lognormal(4.85, 0.4);
  eta_noise  ~ std_normal();
  sigma_x    ~ normal(0, p_sigma_x);
  phiD       ~ normal(0, p_phiD);

  // Likelihood:
  if (inference != 0) {
	for (j in 1:A) target += neg_binomial_2_lpmf( y_data[EpidemicStart:N, j] | E_deathsAge[EpidemicStart:N, j], phiD);
  }// End if
}

generated quantities {
  matrix[N, A] Susceptibles; // Counts of susceptibles at time t, for age group a = 1,..,A.

  vector[N] E_cases;  // Expected infections <lower = 0>
  vector[N] E_deaths; // Expected deaths <lower = 0>

  vector<lower=0>[N] E_infectious;
  vector<lower=0>[N] Rt0;
  
  row_vector[A] convolution_tmp_Rt_Age;
  row_vector[A] Rt_sum_tmp;
  matrix<lower=0>[N, A] effective_EcasesAge;
  vector<lower=0>[N] Rt1;

  matrix[A, N - EpidemicStart + 1] log_like_age; 
  vector[A * (N - EpidemicStart + 1)] log_lik;
  
  real deviance; 
  
  matrix<lower=0>[N, A] Rt_Age = rep_matrix(1.0, N, A);

  //---- Total expected number of new infections per day:
  E_cases = E_casesAge * ones_vector_A;

  //---- Total expected number of new deaths per day:
  E_deaths = E_deathsAge * ones_vector_A;

  for ( i in 1:(N - EpidemicStart + 1) ) {
    for (j in 1:A) log_like_age[j,i] = neg_binomial_2_lpmf(y_data[EpidemicStart + i - 1,j] | E_deathsAge[EpidemicStart + i - 1,j], phiD);
  }// End for
  
  log_lik = matrix_to_vector(log_like_age, A, (N - EpidemicStart + 1));

  //---- Deviance:
  deviance = (-2) * sum(log_lik);
    
  //---- Estimation of R_t:
  E_infectious[1:N0]         = E_cases[1:N0];
  effective_EcasesAge[1:N0,] = E_casesAge[1:N0,];
  Rt1[1:N0]                  = rep_vector(1.0, N0);
  Rt1[(N0+1):N]              = rep_vector(0.0, N-N0);

  for( t in (N0+1):N ) {
	  E_infectious[t] = dot_product( E_cases[1:t-1], tail(GT_rev, t-1) );
	  
	  for (a in 1:A) {
		  
		  for (k in 1:A) convolution_tmp_Rt_Age[k] = (1.0 - (rep_row_vector(1.0, t-1) * E_casesAge[1:(t-1),k] / n_pop[k]) ) * 
												     rho_daily[t,k] * 
												     contact_matrix[a,k];
													 			 
		  Rt_Age[t,a]              = sum( convolution_tmp_Rt_Age );
		  effective_EcasesAge[t,a] = dot_product( E_casesAge[1:t-1,a], tail(GT_rev, t-1) );
		  Rt_sum_tmp[a]            = effective_EcasesAge[t,a] * Rt_Age[t,a];

	  }// End for
	  
	  Rt1[t] = sum(Rt_sum_tmp) / sum(effective_EcasesAge[t,]);
	
  }// End for

  Rt0 = E_cases ./ E_infectious;

  //---- Estimation of age-specific counts of Susceptibles:
  Susceptibles[1,] = E_casesAge[1,];
  for( t in 2:N ) for (a in 1:A) Susceptibles[t,a] = (1.0 - (rep_row_vector(1.0, t-1) * E_casesAge[1:(t-1),a] / n_pop[a]) ) * n_pop[a];
}
