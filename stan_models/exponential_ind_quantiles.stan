




data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0,upper=1>[N] p;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> pv;
  // int<lower=0> n;

}


parameters {
  real<lower=0> lambda;
  real<lower=0> sigma_rho;
}



model {
  lambda ~ normal(0, c);
  sigma_rho ~ normal(0, pv);

  Q ~ normal(-log(1 - p)/lambda, 1/sigma_rho);
}

generated quantities {
  vector[N] pred_q;
  real dist_samp;
  
  for (i in 1:N) {
    pred_q[i] = normal_rng(-log(1 - p[i])/lambda, 1/sigma_rho);
  }
  dist_samp = exponential_rng(lambda);
}

