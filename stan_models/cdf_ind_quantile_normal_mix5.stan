

functions{
      
  
  
  real GM_PDF(vector mu, vector sigma, vector pi, real y) {
    
    real density;
    
    density = pi[1]*exp(normal_lpdf(y | mu[1], sigma[1])) +
              pi[2]*exp(normal_lpdf(y | mu[2], sigma[2])) +
              pi[3]*exp(normal_lpdf(y | mu[3], sigma[3]))// +
              + pi[4]*exp(normal_lpdf(y | mu[4], sigma[4]))
              + pi[5]*exp(normal_lpdf(y | mu[5], sigma[5])) 
	      ;
    
    
    return density;
  }
  
  real GM_CDF(real y, vector mu, vector sigma, vector pi) {
    real prob;
    
    prob = pi[1]*normal_cdf(y | mu[1], sigma[1]) +
           pi[2]*normal_cdf(y | mu[2], sigma[2]) +
           pi[3]*normal_cdf(y | mu[3], sigma[3])
           + pi[4]*normal_cdf(y | mu[4], sigma[4])
           + pi[5]*normal_cdf(y | mu[5], sigma[5])
           ;
           
    return prob;
  }
  
  
}

data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0, upper=1>[N] p;
  int<lower=1> n_components;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd  
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv;
  vector[n_components] alpha; // component weights prior parameter
}



parameters {
  vector[n_components] mus;
  vector<lower=0>[n_components] sigmas;
  real<lower=0> n;
  simplex[n_components] pi;
  real<lower=0> sigma;
  
}

transformed parameters {
  vector<lower=0,upper=1>[N] p_hat;
  
  for (i in 1:N) p_hat[i] = GM_CDF(Q[i], mus, sigmas, pi);
  
}

model {
  pi ~ dirichlet(alpha);
  mus ~ normal(m, c);
  sigmas ~ normal(0, sv);
  sigma ~ normal(0, nv);
  
  
  for (i in 1:N) p_hat[i] ~ normal(p[i], 1/sigma);
}

generated quantities {
  vector[N] pred_q;
  real dist_samp;
  int samp_comp = categorical_rng(pi);
  dist_samp = normal_rng(mus[samp_comp], sigmas[samp_comp]);
}

