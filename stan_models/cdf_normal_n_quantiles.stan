

data {
  int<lower=0> N;
  real<lower=0> n;
  vector[N] Q;
  vector<lower=0,upper=1>[N] p;
  vector[N] inv_Phip;
  cov_matrix[N] QCorr;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
}


parameters {
  real mu;
  real<lower=0> sigma;
}

transformed parameters {
  //vector[N] xi;
  vector<lower=0,upper=1>[N] p_hat;
  for (i in 1:N) {
    //xi[i] = mu + sigma*inv_Phip[i];
    p_hat[i] = normal_cdf(Q[i] | mu, sigma);
  }
  

}



model {
  mu ~ normal(m, c);
  sigma ~ normal(0, sv);

  //Q ~ multi_normal(xi, (sigma^2)*(1/n)*QCorr);
  p_hat ~ multi_normal(p, (1/n)*QCorr);
}

generated quantities {
  //vector[N] pred_q;
  real dist_samp;
  //pred_q = multi_normal_rng(xi, (sigma^2)*(1/n)*QCorr);
  dist_samp = normal_rng(mu, sigma);
  
  vector[N] p_samp;
  vector[N] Q_rep;  // Q such that GM_CDF(Q_rep[i]) ≈ p[i]
  
  p_samp = multi_normal_rng(p, (1/n)*QCorr);
  for (i in 1:N) if (p_samp[i] < 0) p_samp[i] = 1e-6;
  for (i in 1:N) if (p_samp[i] > 1) p_samp[i] = .999999;
  
  for (i in 1:N) {
    Q_rep[i] = mu + sigma * inv_Phi(p_samp[i]); // normal_quantile(p_samp[i], mu, sigma);
  }
}

