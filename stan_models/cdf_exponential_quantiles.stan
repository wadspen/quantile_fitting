
data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0, upper=1>[N] p;
  vector[N] expp;
  cov_matrix[N] QCorr;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}


parameters {
  real<lower=0> lambda;
  real<lower=0> n;
}



transformed parameters {
  vector<lower=0,upper=1>[N] p_hat;
  for (i in 1:N) p_hat[i] = exponential_cdf(Q[i] | lambda);
}


model {
  lambda ~ normal(0, c);
  n ~ normal(0, nv);

  // Q ~ multi_normal(xi, (lambda)*(1/n)*QCorr);
  p_hat ~ multi_normal(p, (1/n)*QCorr);
}

generated quantities {
  // vector[N] pred_q;
  real dist_samp;
  // pred_q = multi_normal_rng(xi, (lambda)*(1/n)*QCorr);
  dist_samp = exponential_rng(lambda);
  
  vector[N] p_samp;
  vector[N] Q_rep;  // Q such that GM_CDF(Q_rep[i]) ≈ p[i]
  
  p_samp = multi_normal_rng(p, (1/n)*QCorr);
  for (i in 1:N) if (p_samp[i] < 0) p_samp[i] = 1e-6;
  for (i in 1:N) if (p_samp[i] > 1) p_samp[i] = .999999;
  
  for (i in 1:N) {
    Q_rep[i] = -log(1 - p_samp[i])/lambda; // normal_quantile(p_samp[i], mu, sigma);
  }
}

