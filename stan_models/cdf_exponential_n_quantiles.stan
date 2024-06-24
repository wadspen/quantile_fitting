
data {
  int<lower=0> N;
  real<lower=0> n;
  vector<lower=0, upper=1>[N] p;
  vector[N] Q;
  cov_matrix[N] QCorr;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}


parameters {
  real<lower=0> lambda;
}

transformed parameters {
  vector<lower=0,upper=1>[N] p_hat;
  for (i in 1:N) p_hat[i] = exponential_cdf(Q[i] | lambda);
}


model {
  lambda ~ normal(0, c);

  p_hat ~ multi_normal(p, (1/n)*QCorr);
}

generated quantities {
  // vector[N] pred_q;
  real dist_samp;
  // pred_q = multi_normal_rng(xi, (1/lambda)*(1/n)*QCorr);
  dist_samp = exponential_rng(lambda);
}

