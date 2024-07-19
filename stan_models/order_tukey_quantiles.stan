


functions {
  real orderstatistics (real n, int N, vector p, vector U){
    real lpdf = 0;
      lpdf += lgamma(n + 1) - lgamma(n*p[1]) - lgamma(n - n*p[N] + 1) ;
      lpdf += (n*p[1] - 1)*log(U[1]);
      lpdf += (n - n*p[N])*log(1 - U[N]) ;
      for (i in 2:N){
        lpdf += -lgamma (n*p[i] - n*p[i - 1]) ;
        lpdf += (n*p[i]-n*p[i - 1] - 1) * log(U[i] - U[i - 1]) ;
      }
      return lpdf ;
    }
    
    real tukey_q(real logit_p, real lambda) {
        real log_p = log_inv_logit(logit_p);
        real log1m_p = log1m_inv_logit(logit_p);
        real s = lambda*log_p;
        real t = lambda*log1m_p;
        if (abs(s) > 1e-6 || abs(t) > 1e-6)
          return (exp(s) - exp(t))/lambda;
        else
          return (log_p - log1m_p) + (s*log_p - t*log1m_p);
    }
    real tukey_rng(real lambda) {
        real q = uniform_rng(0, 1);
        return tukey_q(logit(q), lambda);
    }
    vector tukey_system(vector logit_p, vector theta, real[] y, int[] z) {
        return [tukey_q(logit_p[1], theta[1])-theta[2]]';
    }
    real tukey_lcdf(real x, real lambda) {
        vector[1] logit_p = algebra_solver(tukey_system, [0.0]', [lambda, x]', {0.0},{0});
        return log_inv_logit(logit_p[1]);
    }
    real tukey_lpdf(real x, real lambda) {
        vector[1] logit_p = algebra_solver(tukey_system, [0.0]', [lambda, x]', {0.0},{0});
        real log_p = log_inv_logit(logit_p[1]);
        real log1m_p = log1m_inv_logit(logit_p[1]);
        return -log_sum_exp((lambda-1)*log_p, (lambda-1)*log1m_p);
    }
}

data {
  int N; // number of observed quantiles
  vector[N] Q; // quantiles
  vector[N] p; // probabilities
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}

parameters {
  real lambda;
  real<lower=0> n;
}

transformed parameters {
  vector[N] U;
  for (i in 1:N)
    U[i] = exp(tukey_lcdf(Q[i] | lambda));
}

model {
  lambda ~ normal(0, c);
  n ~ normal(0, nv);

  
  target += orderstatistics(n , N, p, U);
  
  for (i in 1:N)
    target += tukey_lpdf (Q[i] | lambda);
}

generated quantities {
  real dist_samp = tukey_rng(lambda);
  // real log_prob = orderstatistics(n, N, p, U);
  // for (i in 1:N)
  //   log_prob += normal_lpdf(Q[i] | mu, sigma);
}

