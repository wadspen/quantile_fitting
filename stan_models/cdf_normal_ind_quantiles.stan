




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
  
  vector[N] p_samp;
  vector[N] Q_rep;  // Q such that GM_CDF(Q_rep[i]) ≈ p[i]
  
  for (i in 1:N) {
    p_samp[i] = normal_rng(p[i], 1/sigma_rho);
    if (p_samp[i] < 0) p_samp[i] = 1e-6;
    if (p_samp[i] > 1) p_samp[i] = .999999;
  }
  
  for (i in 1:N) {
    Q_rep[i] = mu + sigma * inv_Phi(p_samp[i]); // normal_quantile(p_samp[i], mu, sigma);
  }
}

