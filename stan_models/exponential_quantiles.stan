
data {
  int<lower=0> N;
  vector[N] Q;
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
  vector[N] xi;
  
  for (i in 1:N) {
    xi[i] = expp[i]/lambda;
  }

}


model {
  lambda ~ normal(0, c);
  n ~ normal(0, nv);

  Q ~ multi_normal(xi, (1/lambda)*(1/n)*QCorr);
}

generated quantities {
  vector[N] pred_q;
  real dist_samp;
  pred_q = multi_normal_rng(xi, (1/lambda)*(1/n)*QCorr);
  dist_samp = exponential_rng(lambda);
}

