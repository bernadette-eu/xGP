functions {
  
  matrix compute_covariance_matrix(int M, int[] start_idx, int[] end_idx,
                                   matrix min_dist_times, 
                                   real[] sigma_A,
                                   real sigma_M) {
    int N = rows(min_dist_times);
    matrix[N, N] K = rep_matrix(0.0, N, N); // final covariance matrix

    // Add shared M-level covariance to entire matrix
    for (i in 1:N) {
      for (j in i:N) {
        real v = sigma_M^2 * min_dist_times[i, j];
        K[i, j] = v;
        if (i != j) K[j, i] = v; // enforce symmetry
      }
    }

    // Add group-specific (diagonal block) covariance
    for (k in 1:M) {
      real sigma2_A = square(sigma_A[k]);
      for (i in start_idx[k]:end_idx[k]) {
        for (j in i:end_idx[k]) {
          real v = sigma2_A * min_dist_times[i, j];
          K[i, j] += v;
          if (i != j) K[j, i] += v;
        }
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
  real<lower=0> length_scale_A[M];
  real<lower=0> sigma_A[M];
  real<lower=0> length_scale_M;
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

