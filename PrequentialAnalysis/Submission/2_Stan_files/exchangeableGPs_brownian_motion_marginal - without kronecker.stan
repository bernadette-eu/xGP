functions {
  
  matrix compute_covariance_matrix(int M, int[] start_idx, int[] end_idx, 
                                  matrix min_dist_times, 
                                  real sigma_A, 
                                  real sigma_M) {
    int N = rows(min_dist_times);  // Total number of observations
    matrix[N, N] K = rep_matrix(0.0, N, N);  // Initialize covariance matrix
    
    // Precompute constants
    real sigma_sq_A = square(sigma_A);
    real sigma_sq_M = square(sigma_M);
    
    // Fill diagonal blocks
    for (k in 1:M) {
        int start_k = start_idx[k];
        int end_k = end_idx[k];
        int n_k = end_k - start_k + 1;
        
        // Extract diagonal block
        matrix[n_k, n_k] sq_block = min_dist_times[start_k:end_k, start_k:end_k];
        
        // Covariance: Sum of kernels for diagonal blocks
        matrix[n_k, n_k] K_block = sigma_sq_A * sq_block +
                                   sigma_sq_M * sq_block;
        
        // Assign back to main matrix
        K[start_k:end_k, start_k:end_k] = K_block;
    }
    
    // Fill off-diagonal blocks
    for (k in 1:(M - 1)) {
        for (m in (k + 1):M) {
            int start_k = start_idx[k];
            int end_k = end_idx[k];
            int start_m = start_idx[m];
            int end_m = end_idx[m];
            
            // Extract off-diagonal block
            matrix[end_k - start_k + 1, end_m - start_m + 1] 
              sq_block = min_dist_times[start_k:end_k, start_m:end_m];
            
            // Covariance: Single kernel for off-diagonal blocks
            matrix[end_k - start_k + 1, end_m - start_m + 1] 
              K_block = sigma_sq_M  * sq_block;
            
            // Assign back to main matrix (both symmetric parts)
            K[start_k:end_k, start_m:end_m] = K_block;
            K[start_m:end_m, start_k:end_k] = K_block';  // Mirror
        }
    }
    
    return K;
}

}
data{
  int<lower=0> M; //number of groups-countries-individuals
  int<lower=0>  dim_Y[M]; //number of timepoints-obervations per group-country-individual
  int<lower=0> N; //Total number of points
  
  vector[N] Y; //Concatenated data points, we get Y as the vector of all the observed directories stacked one after the other 

  vector<lower=0>[N] Times; //Concatenated time points, we get Times as the vector of all the observed times stacked one after the other 


  
}
transformed data {
  int<lower=0>  start_idx[M]; //Vector indicating where its index starts for each group-country
  int<lower=0> end_idx[M];//Vector indicating where its index ends for each group-country

  vector[N] mu = rep_vector(0, N);
  
  matrix[N,N] I_n=diag_matrix(rep_vector(1.0, N));
  real delta = 1e-9;
  matrix<lower=0>[N,N] min_dist_times;
  
  start_idx[1] = 1;
  end_idx[1] = dim_Y[1];
  for (i in 2:M) {
    start_idx[i] = end_idx[i - 1] + 1;
    end_idx[i] = end_idx[i - 1] + dim_Y[i];
  }
  

  for (i in 1:N) {
    for(j in 1:N){
        min_dist_times[i,j] = fmin(Times[i],Times[j]);

    }
  }

}
parameters {
  real<lower=0> sigma_A;
  real<lower=0> sigma_M; 
  real<lower=0> sigma_e;
}
transformed parameters{
  matrix[N, N] L_K;
  matrix[N, N] K;


  K=compute_covariance_matrix(M,start_idx,end_idx,min_dist_times,sigma_A,sigma_M);

  L_K= cholesky_decompose(K+I_n*(square(sigma_e)+delta));
}
model {
  //priors
  sigma_A~std_normal();
  sigma_M~std_normal(); 
  sigma_e~std_normal();
  //marginal likelihood
  Y ~ multi_normal_cholesky(mu, L_K);
}


