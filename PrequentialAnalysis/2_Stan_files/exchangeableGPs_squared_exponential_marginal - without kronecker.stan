functions {
  // matrix compute_covariance_matrix(int M, int[] start_idx, int[] end_idx, matrix sq_dist_times, 
  //                           real sigma_A, real length_scale_A, real sigma_M, real length_scale_M) {
  //   int N = rows(sq_dist_times); // Total number of observations
  //   matrix[N, N] K; // Initialize the covariance matrix
  // 
  //   // Fill the covariance matrix
  //   for (k in 1:M) {
  //     for (m in 1:M) {
  //       for (i in start_idx[k]:end_idx[k]) {
  //         for (j in start_idx[m]:end_idx[m]) {
  //           if (k == m) { // Diagonal blocks
  //             K[i, j] = sigma_A^2 * exp(-0.5 * sq_dist_times[i,j] / length_scale_A^2) +
  //                        sigma_M^2 * exp(-0.5 * sq_dist_times[i,j] / length_scale_M^2);
  //           } else { // Off-diagonal blocks
  //             K[i, j] = sigma_M^2 * exp(-0.5 * sq_dist_times[i,j] / length_scale_M^2);
  //           }
  //         }
  //       }
  //     }
  //   }
  // 
  //   return K;
  // }
  
  matrix compute_covariance_matrix(int M, int[] start_idx, int[] end_idx, 
                                  matrix sq_dist_times, 
                                  real sigma_A, real length_scale_A, 
                                  real sigma_M, real length_scale_M) {
    int N = rows(sq_dist_times);  // Total number of observations
    matrix[N, N] K = rep_matrix(0.0, N, N);  // Initialize covariance matrix
    
    // Precompute constants
    real inv_len_sq_A = -0.5 / square(length_scale_A);
    real inv_len_sq_M = -0.5 / square(length_scale_M);
    real sigma_sq_A = square(sigma_A);
    real sigma_sq_M = square(sigma_M);
    
    // Fill diagonal blocks
    for (k in 1:M) {
        int start_k = start_idx[k];
        int end_k = end_idx[k];
        int n_k = end_k - start_k + 1;
        
        // Extract diagonal block
        matrix[n_k, n_k] sq_block = sq_dist_times[start_k:end_k, start_k:end_k];
        
        // Covariance: Sum of kernels for diagonal blocks
        matrix[n_k, n_k] K_block = sigma_sq_A * exp(inv_len_sq_A * sq_block) +
                                   sigma_sq_M * exp(inv_len_sq_M * sq_block);
        
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
              sq_block = sq_dist_times[start_k:end_k, start_m:end_m];
            
            // Covariance: Single kernel for off-diagonal blocks
            matrix[end_k - start_k + 1, end_m - start_m + 1] 
              K_block = sigma_sq_M * exp(inv_len_sq_M * sq_block);
            
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
  matrix<lower=0>[N,N] sq_dist_times;
  
  start_idx[1] = 1;
  end_idx[1] = dim_Y[1];
  for (i in 2:M) {
    start_idx[i] = end_idx[i - 1] + 1;
    end_idx[i] = end_idx[i - 1] + dim_Y[i];
  }
  

  for (i in 1:N) {
    sq_dist_times[,i ] = square(Times[i] - Times);
  }

}
parameters {
  real<lower=0> length_scale_A;
  real<lower=0> sigma_A;
  real<lower=0> length_scale_M;
  real<lower=0> sigma_M; 
  real<lower=0> sigma_e;
}
transformed parameters{
  matrix[N, N] L_K;
  matrix[N, N] K;


  K=compute_covariance_matrix(M,start_idx,end_idx,sq_dist_times,sigma_A,length_scale_A,sigma_M,length_scale_M);

  L_K= cholesky_decompose(K+I_n*(square(sigma_e)+delta));
}
model {
  //priors
  length_scale_A~std_normal();
  sigma_A~std_normal();
  length_scale_M~std_normal();
  sigma_M~std_normal(); 
  sigma_e~std_normal();
  //marginal likelihood
  Y ~ multi_normal_cholesky(mu, L_K);
}

//generated quantities{
//posterior predictive
//vector[N_pred] Y_pred;
//vector[N_pred] mu_Y_pred;
//matrix[N_pred, N_pred] K_Y_pred;
//matrix[num_T_pred,num_T] K_A_test;
//matrix[num_T_pred,num_T] K_M_test;
//matrix[num_T_pred,num_T_pred] K_A_test_self;
//matrix[num_T_pred,num_T_pred] K_M_test_self;
//matrix[N_pred,N_pred] K_test_self;
//matrix[N_pred,N] K_test;

//f
//vector[N] f;
//vector[N] mu_f_given_Y;
//matrix[N, N] K_f_given_Y;


//fitted
//vector[N] Y_fitted;

//f_M
  //vector[num_T] f_M;
  //matrix[num_T,M] matrix_f;
  //vector[num_T] sum_f_per_M;
  //vector[num_T] mu_f_M_given_f;
  //matrix[num_T, num_T] K_f_M_given_f;
  
  //matrix[N, N] inv_K;
  //matrix[num_T, num_T] inv_K_A;
//matrix[num_T, num_T] inv_K_M;

//log_likelihood for WAIC
//real log_lik[N];

 



//inv_K=inverse_spd(K);
//inv_K_A=inverse_spd(K_A);
//inv_K_M=inverse_spd(K_M);
//f
//mu_f_given_Y=transpose(K-I_n*square(sigma_e))*inv_K*Y;
//K_f_given_Y=(K-I_n*square(sigma_e))-transpose(K-I_n*square(sigma_e))*inv_K*(K-I_n*square(sigma_e));
//f=multi_normal_rng(mu_f_given_Y,K_f_given_Y);
//fitted values y hat
//Y_fitted=mu_f_given_Y;
//f_M
//K_f_M_given_f=inverse_spd(inv_K_M+M*inv_K_A);
//matrix_f= to_matrix(f, num_T, M);
//  for (i in 1:num_T) {
//    sum_f_per_M[i] = sum(matrix_f[i,1:M]);  // Sum each row
//  }

//mu_f_M_given_f=K_f_M_given_f*inv_K_A*sum_f_per_M;
//f_M=multi_normal_rng(mu_f_M_given_f,K_f_M_given_f);

//posterior predictive
//K_A_test=gp_exp_quad_cov(Times_pred,Times, sigma_A, length_scale_A);
//K_M_test=gp_exp_quad_cov(Times_pred,Times, sigma_M, length_scale_M);
//K_A_test_self=gp_exp_quad_cov(Times_pred, sigma_A, length_scale_A);
//K_M_test_self=gp_exp_quad_cov(Times_pred, sigma_M, length_scale_M);
 
///K_test=kronecker_product(Ones,K_M_test)+kronecker_product(I_m,K_A_test);
//K_test_self=kronecker_product(Ones,K_M_test_self)+kronecker_product(I_m,K_A_test_self);
//mu_Y_pred=K_test*inv_K*Y;
//K_Y_pred=K_test_self-K_test*inv_K*transpose(K_test);

//print(K_A);
//print(dims(K_A));

//Y_pred=multi_normal_rng(mu_Y_pred,K_Y_pred);


//loglikelihood
//for(i in 1:N){
//  log_lik[i]=normal_lpdf(Y[i] | f[i], sigma_e);
//}
//
//}
