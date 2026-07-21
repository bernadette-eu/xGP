functions {
  vector matrix_to_vector(matrix X, int num_rows, int num_cols) {
    return to_vector(to_row_vector(X)');
  }
}

data {
	int<lower = 1> A;                  
	int<lower = 1> N;                  
	int EpidemicStart;
	int y_data[N, A];                 
}

parameters {
	real<lower = 0> phiD;              
	matrix<lower = 0>[N, A] E_deathsAge;	
}

model {
}

generated quantities {
  matrix[N - EpidemicStart + 1, A] log_like_age; 
  vector[A * (N - EpidemicStart + 1)] log_lik;
  real deviance;

  for ( i in 1:(N - EpidemicStart + 1) ) {
	for (j in 1:A) log_like_age[i,j] = neg_binomial_2_lpmf(y_data[EpidemicStart + i - 1,j] | E_deathsAge[EpidemicStart + i - 1,j], phiD);
  }// End for

  log_lik = matrix_to_vector(log_like_age, A, (N - EpidemicStart + 1));

  //---- Deviance:
  deviance = (-2) * sum(log_lik);
}
