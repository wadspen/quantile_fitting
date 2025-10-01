


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
  int N; // number of observed quantiles
  real<lower=0> n;
  vector[N] Q; // quantiles
  vector[N] p; // probabilities
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}

parameters {
  real mu;
  real<lower=0> sigma;

}

transformed parameters {
  vector[N] U;
  for (i in 1:N)
    U[i] = normal_cdf(Q[i] | mu, sigma);
}

model {
  mu ~ normal(m, c);
  sigma ~ normal(0, sv);

  
  target += orderstatistics(n , N, p, U);
  
  for (i in 1:N)
    target += normal_lpdf (Q[i] | mu, sigma);
}

generated quantities {
  real dist_samp = normal_rng(mu, sigma);
  // real log_prob = orderstatistics(n, N, p, U);
  // for (i in 1:N)
  //   log_prob += normal_lpdf(Q[i] | mu, sigma);
}

