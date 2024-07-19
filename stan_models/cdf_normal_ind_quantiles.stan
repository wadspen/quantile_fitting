




data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0,upper=1>[N] p;
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

transformed parameters {
  vector<lower=0,upper=1>[N] p_hat;
  for (i in 1:N) p_hat[i] = normal_cdf(Q[i] | mu, sigma);
}



model {
  mu ~ normal(m, c);
  sigma ~ normal(0, sv); 
  sigma_rho ~ normal(0, pv);
  
  p_hat ~ normal(p, 1/sigma_rho);
  //Q ~ normal(mu + sigma*inv_Phip, 1/sigma_rho);
}

generated quantities {
  //vector[N] pred_q;
  real dist_samp;
  
  //for (i in 1:N) {
  //  pred_q[i] = normal_rng(mu + sigma*inv_Phip[i], sigma_rho);
  //}
  dist_samp = normal_rng(mu, sigma);
}

