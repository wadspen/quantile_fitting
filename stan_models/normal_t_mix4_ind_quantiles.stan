

functions{
      
  
  
  real GM_PDF(vector mu, vector sigma, vector pi, real y) {
    
    real density;
    
    density = pi[1]*exp(normal_lpdf(y | mu[1], sigma[1])) +
              pi[2]*exp(normal_lpdf(y | mu[2], sigma[2])) +
              pi[3]*exp(normal_lpdf(y | mu[3], sigma[3]))// +
              + pi[4]*exp(normal_lpdf(y | mu[4], sigma[4]))
              ;
    
    
    return density;
  }
  
  real GM_CDF(real y, vector mu, vector sigma, vector pi) {
    real prob;
    
    prob = pi[1]*normal_cdf(y | mu[1], sigma[1]) +
           pi[2]*normal_cdf(y | mu[2], sigma[2]) +
           pi[3]*normal_cdf(y | mu[3], sigma[3])
           + pi[4]*normal_cdf(y | mu[4], sigma[4])
           ;
           
    return prob;
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
  vector[n_components] alpha;
}



parameters {
  vector[n_components] mus;
  vector<lower=0>[n_components] sigmas;
  real<lower=0> n;
  simplex[n_components] pi;
  real<lower=0> sigma;
  
  
  
}

transformed parameters {
  // vector[n_components] pit;
  vector<lower=0,upper=1>[N] p_hat;
  
  for (i in 1:N) p_hat[i] = GM_CDF(Q[i], mus, sigmas, pi);
  

  // pit = softmax(pi);
  


}

model {
  // pi ~ normal(4,5);
  pi ~ dirichlet(alpha);
  mus ~ normal(4, 5);
  sigmas ~ normal(0,7);
  sigma ~ normal(0,1);
  
  // for (i in 1:N) Q[i] ~ normal(Qi[i], sigma);
  for (i in 1:N) p_hat[i] ~ normal(p[i], sigma);
  //Q ~ multi_normal(Qi, diag_matrix(square(sigma)));
}

generated quantities {
  vector[N] pred_q;
  real dist_samps;
  // for (i in 1:N) pred_q[i] = normal_rng(Qi[i], sigma);
  //pred_q = multi_normal_rng(Qi, diag_matrix(square(sigma)));
  int samp_comp = categorical_rng(pi);
  dist_samps = normal_rng(mus[samp_comp], sigmas[samp_comp]);
}

