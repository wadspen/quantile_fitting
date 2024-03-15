




functions {
  real normal_quantile(real p, real mu, real sigma) {
    if (is_nan(p) || is_inf(p) || p < 0 || p > 1)
    reject("normal_quantile: p must be finite and between 0 and 1; ", 
    "found p = ", p);
    
    real q;
    
    q = mu + inv_Phi(p)*(sigma);
    return q;
    
  }
  
  real lognormal_quantile(real p, real mu, real sigma) {
    if (is_nan(p) || is_inf(p) || p < 0 || p > 1)
    reject("normal_quantile: p must be finite and between 0 and 1; ", 
    "found p = ", p);
    
    real q;
    
    q = exp(mu + inv_Phi(p)*(sigma));
    return q;
    
  }

}


data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0, upper=1>[N] p;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
  // int<lower=0> n;
}


parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> n;
}

transformed parameters {
  vector[N] xi;
  matrix<lower=0>[N, N] q_var;

  for (i in 1:N) {
    xi[i] = normal_quantile(p[i], mu, sigma);
  }
  
  for (i in 1:N) {
    for (j in 1:i) {
      
      q_var[i, j] = (2*pi()*(sigma^2)*fmin(p[i],p[j])*(1 - fmax(p[i],p[j]))) / 
                    (n*exp((-.5)*(inv_Phi(p[i])^2 + inv_Phi(p[j])^2)));
      
      q_var[j, i] = q_var[i, j];
      
    }
  }

}


model {
  mu ~ normal(m, c);
  sigma ~ normal(0, sv);
  n ~ normal(0, nv);

  Q ~ multi_normal(xi, q_var);
}

generated quantities {
  vector[N] pred_q;
  real dist_samps;
  pred_q = multi_normal_rng(xi, q_var);
  dist_samps = normal_rng(mu, sigma);
}

