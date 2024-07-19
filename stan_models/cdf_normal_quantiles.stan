
data {
  int<lower=0> N;
  vector[N] Q;
  vector[N] inv_Phip;
  vector<lower=0,upper=1>[N] p;
  cov_matrix[N] QCorr;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}


parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> n;
}

transformed parameters {
  vector<lower=0,upper=1>[N] p_hat;
  
  for (i in 1:N) {
    p_hat[i] = normal_cdf(Q[i] | mu, sigma);
  }

}


model {
  mu ~ normal(m, c);
  sigma ~ normal(0, sv);
  n ~ normal(0, nv);

  //Q ~ multi_normal(xi, (sigma^2)*(1/n)*QCorr);
  p_hat ~ multi_normal(p, (1/n)*QCorr);
}

generated quantities {
  //vector[N] pred_q;
  real dist_samp;
  //pred_q = multi_normal_rng(xi, (sigma^2)*(1/n)*QCorr);
  dist_samp = normal_rng(mu, sigma);
}

