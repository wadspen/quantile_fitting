//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//





// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector<lower=0, upper=1>[N] p;
  ordered[N] Q;
  vector[N] inv_Phip;
  cov_matrix[N] QCorr;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real a1;
  real a2;
  real a3;
  real a4;
  real a5;
  real a6;
  real a7;
  real a8;
  // real a9;
  // real a10;
  real<lower=0> n;
  
}

transformed parameters {
  vector[N] mu;
  vector[N] log_sigma;
  cov_matrix[N] QRho;
  
  mu = a1 + a4*(p - .5)  
  + a5*(p - .5)^2
  + a7*(p - .5)^3
  // + a9*(p - .5)^4
  ;
  log_sigma = a2 + a3*(p - .5) 
  + a6*(p - .5)^2
  + a8*(p - .5)^3
  // + a10*(p - .5)^4
  ;
  
  QRho = quad_form_diag(QCorr, exp(log_sigma))*(n^(-1));
  

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a1 ~ normal(0, 30);
  a2 ~ normal(0, 30);
  a3 ~ normal(0, 30);
  a4 ~ normal(0, 30);
  a5 ~ normal(0, 30);
  a6 ~ normal(0, 30);
  a7 ~ normal(0, 30);
  a8 ~ normal(0, 30);
  // a9 ~ normal(0, 30);
  // a10 ~ normal(0, 30);
  n ~ normal(0, 5000);
  
  Q ~ multi_normal(mu + exp(log_sigma).*inv_Phip, QRho);
  
}

generated quantities {
  vector[N] pred_q;
  real log_score;
  real dist_samp;
  
  pred_q = multi_normal_rng(mu + exp(log_sigma).*inv_Phip, (n^(-1))*QRho);
  
  real<lower=0,upper=1> samp_p = uniform_rng(0,1); 
  
  dist_samp = a1 + a4*(samp_p - .5)  
                + a5*(samp_p - .5)^2
                + a7*(samp_p - .5)^3
                // + a9*(samp_p - .5)^4
                + exp(a2 + a3*(samp_p - .5) 
                + a6*(samp_p - .5)^2
                + a8*(samp_p - .5)^3
                // + a10*(samp_p - .5)^4
                ).*inv_Phi(samp_p);
}

