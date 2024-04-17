




data {
  int<lower=0> N;
  vector[N] Q;
  vector[N] inv_Phip;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> pv;
  // int<lower=0> n;
}


parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> sigma_rho;
}



model {
  mu ~ normal(m, c);
  sigma ~ normal(0, sv);
  sigma_rho ~ normal(0, pv);

  Q ~ normal(mu + sigma*inv_Phip, sigma_rho);
}

generated quantities {
  vector[N] pred_q;
  real dist_samps;
  
  for (i in 1:N) {
    pred_q[i] = normal_rng(mu + sigma*inv_Phip[i], sigma_rho);
  }
  dist_samps = normal_rng(mu, sigma);
}

