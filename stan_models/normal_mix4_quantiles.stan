

functions{
      
  vector GMInv_CDF(vector Q, data real p, vector pi, vector mu, 
                   vector sigma, int N) {
    vector[1] z;
    real cdf;
    // this encodes the cdf from which we want to sample
    cdf = pi[1]*exp(normal_lcdf(Q | mu[1], sigma[1])) +
          pi[2]*exp(normal_lcdf(Q | mu[2], sigma[2])) +
          pi[3]*exp(normal_lcdf(Q | mu[3], sigma[3])) //+
          //pi[4]*exp(normal_lcdf(Q | mu[4], sigma[4]))
          ;
    z[1] = cdf - p;
    return z;
  }
  
  
  real GM_PDF(vector mu, vector sigma, vector pi, real y) {
    
    real density;
    
    density = pi[1]*(exp(normal_lpdf(y | mu[1], sigma[1]))) +
              pi[2]*(exp(normal_lpdf(y | mu[2], sigma[2]))) +
              pi[3]*(exp(normal_lpdf(y | mu[3], sigma[3]))) +
              //pi[4]*exp(normal_lpdf(y | mu[4], sigma[4])) +
              .00001;
              
    // density = log_sum_exp({log(pi[1]) + normal_lpdf(y | mu[1], sigma[1]),
    //                        log(pi[2]) + normal_lpdf(y | mu[2], sigma[2]),
    //                        log(pi[3]) + normal_lpdf(y | mu[3], sigma[3]),
    //                        log(pi[4]) + normal_lpdf(y | mu[4], sigma[4]),
    //                        log(pi[5]) + normal_lpdf(y | mu[5], sigma[5])});
    
    return density;
  }
  
  
}

data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0, upper=1>[N] p;
  int<lower=1> n_components;
  real m;
  real<lower=0> c; 
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}

transformed data{
  real scaling_step = 1e-5;
  real rel_tol = 1e-6;
  real f_tol = 1;
  int max_steps = 3000;
  vector[1] y_guess = [0.5]';
}

parameters {
  ordered[n_components] mus;
  vector<lower=0>[n_components] sigmas;
  real<lower=0> n;
  vector<lower=0,upper=1>[n_components] pi;

}

transformed parameters {
  vector[N] Qi;
  cov_matrix[N] q_var;
  vector[n_components] pit;
  

  pit = softmax(pi);
  
  for (i in 1:N) {
    // // vector[1] Q_guess = [Q[i]]';
    Qi[i] = solve_powell_tol(GMInv_CDF, y_guess,
                             rel_tol, f_tol, max_steps,
                             p[i], pit, mus, sigmas, N)[1];
                             
    // Qi[i] = solve_newton_tol(GMInv_CDF, y_guess,
    //                          scaling_step, f_tol, max_steps,
    //                          p[i], pit, mus, sigmas, N)[1];
  }
  
  for (i in 1:N) {
    for (j in 1:i) {
      
      q_var[i, j] = ((fmin(p[i],p[j])*(1 - fmax(p[i],p[j]))) / 
                    (n*GM_PDF(mus, sigmas, pit, Qi[i])*
                     GM_PDF(mus, sigmas, pit, Qi[j]) + .000001)) +
                     .000001;
      
      q_var[j, i] = q_var[i, j];
      
    }
  }
  
  // q_var_cholesky = cholesky_decompose(q_var);

}

model {
  pi ~ normal(4,5);
  mus ~ normal(4, 5);
  sigmas ~ normal(0,7);
  n ~ normal(0, nv);

  Q ~ multi_normal(Qi, q_var);
  // Q ~ multi_normal_cholesky(Qi, q_var_cholesky);
}

generated quantities {
  vector[N] pred_q;
  // real log_score;
  real dist_samps;
  // for (i in 1:N) cdf[i] = normal_rng(Q[i], mu, sigma);
  pred_q = multi_normal_rng(Qi, q_var);
  // dist_samps = normal_rng(mu, sigma);
  // dist_samps = unif_rng(0,1)
  // real unif_samp = unif_rng(0,1);
  int samp_comp = categorical_rng(pit);
  dist_samps = normal_rng(mus[samp_comp], sigmas[samp_comp]);
}

