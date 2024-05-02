


functions {
  real orderstatistics (real n, int N, vector p, vector U){
    real lpdf = 0;
      lpdf += lgamma(n + 1) - lgamma(n*p[1]) - lgamma(n - n*p[N] + 1) ;
      lpdf += (n*p[1] - 1)*log(U[1]);
      lpdf += (n - n*p[N])*log(1 - U[N]) ;
      for (i in 2:N){
        lpdf += -lgamma (n*p[i] - n*p[i - 1]) ;
        lpdf += (n*p[i]-n*p[i - 1] - 1) * log(U[i] - U[i - 1]) ;
      }
      return lpdf ;
    }
}

data {
  int<lower=0> N;
  int<lower=0> n;
  vector[N] Q;
  vector<lower=0, upper=1>[N] p;
  int<lower=1> n_components;
  real m;
  real<lower=0> c; 
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}

parameters {
  ordered[n_components] mus;
  vector<lower=0>[n_components] sigmas;
  vector<lower=0,upper=1>[n_components] pi;
}

transformed parameters {
  vector[N] U;
  vector[n_components] pit;
  
  pit = softmax(pi);
  
  for (i in 1:N)
    U[i] = pit[1]*normal_cdf(Q[i]| mus[1], sigmas[1]) + 
           pit[2]*normal_cdf(Q[i]| mus[2], sigmas[2]) +
           pit[3]*normal_cdf(Q[i]| mus[3], sigmas[3]) //+
           //pit[4]*normal_cdf(Q[i]| mus[4], sigmas[4]);
           ;
}

model {
  pi ~ normal(4,5);
  mus ~ normal(4, 5);
  sigmas ~ normal(0,7);

  
  target += orderstatistics(n, N, p, U);
  
  for (i in 1:N)
    target += log_sum_exp({log(pit[1]) + normal_lpdf(Q[i] | mus[1], sigmas[1]),
                          log(pit[2]) + normal_lpdf(Q[i] | mus[2], sigmas[2]),
                          log(pit[3]) + normal_lpdf(Q[i] | mus[3], sigmas[3])});//,
                          //log(pit[4]) + normal_lpdf(Q[i] | mus[4], sigmas[4])});
    
}

generated quantities {
  vector[N] pred_q;
  real dist_samps;
  int samp_comp = categorical_rng(pit);
  dist_samps = normal_rng(mus[samp_comp], sigmas[samp_comp]);
}

