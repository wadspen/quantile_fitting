functions {
    real tukey_q(real p, real lambda) {
        // real log_p = log_inv_logit(logit_p);
        // real log1m_p = log1m_inv_logit(logit_p);
        // real s = lambda*log_p;
        // real t = lambda*log1m_p;
        // if (abs(s) > 1e-6 || abs(t) > 1e-6)
        //   return (exp(s) - exp(t))/lambda;
        // else
        //   return (log_p - log1m_p) + (s*log_p - t*log1m_p);
        
        real x = (1/lambda) * (p^lambda - (1 - p)^lambda);
        return x;
    }
    real tukey_qd(real p, real lambda) {
        // real log_p = log_inv_logit(logit_p);
        // real log1m_p = log1m_inv_logit(logit_p);
        // real s = (lambda-1)*log_p;
        // real t = (lambda-1)*log1m_p;
        // if (abs(s) > 1e-6 || abs(t) > 1e-6)
        //   return (exp(s) - exp(t));
        // else
        //   return (log_p - log1m_p) + (s*log_p - t*log1m_p);
        
        real x = (p^(lambda - 1) + (1 - p)^(lambda - 1));
        return x;
    }
    real tukey_rng(real lambda) {
        real q = uniform_rng(0, 1);
        return tukey_q(q, lambda);
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
  int<lower=0> N;
  vector[N] Q;
  vector[N] p;
  cov_matrix[N] QCorr;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
}


parameters {
  real lambda;
  real<lower=0> n;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] xi;
  cov_matrix[N] q_var;
  
  for (i in 1:N) {
    xi[i] = tukey_q(p[i], lambda);
    // xi[i] = (1/lambda) * (p[i]^lambda - (1 - p[i])^lambda);
  }
  
  
  for (i in 1:N) {
    for (j in 1:i) {

      q_var[i, j] = (QCorr[i,j] *
                    (tukey_qd(p[i], lambda) *
                      tukey_qd(p[j], lambda)));

      q_var[j, i] = q_var[i, j];

      if (i == j) {q_var[i,j] = q_var[i,j] + .0000001;}

    }
  }

}


model {
  lambda ~ normal(0, c);
  n ~ normal(0, nv);
  sigma ~ normal(0, 7);
  // for (i in 1:N) {Q[i] ~ normal(xi[i], sigma);}
  Q ~ multi_normal(xi, (1/n)*q_var);
}

generated quantities {
  vector[N] pred_q;
  real dist_samp;
  // for (i in 1:N) {pred_q[i] = normal_rng(xi[i], sigma);}
  pred_q = multi_normal_rng(xi, (1/n)*q_var);
  dist_samp = tukey_rng(lambda);
}

