

functions{
      
  vector GMInv_CDF(vector Q, data real p, vector pi, vector mu, 
                   vector sigma, int N, int num_comp) {
    vector[1] z;
    real cdf;
    vector[num_comp] comp_cdf;
    for (n in 1:num_comp) {
      comp_cdf[n] = pi[n]*exp(normal_lcdf(Q | mu[n], sigma[n]));
    }
    cdf = sum(comp_cdf);
    z[1] = cdf - p;
    return z;
  }
  
  real GM_CDF(real y, vector mu, vector sigma, vector pi, int num_comp) {
    real prob;
    vector[num_comp] comp_prob;
    
    for (n in 1:num_comp) {
      comp_prob[n] = pi[n]*normal_cdf(y | mu[n], sigma[n]);
    }
    
    prob = sum(comp_prob);
           
    return prob;
  }
  
  real GM_PDF(vector mu, vector sigma, vector pi, real y, int num_comp) {
    
    real density;
    vector[num_comp] comp_dens;
    
    for (n in 1:num_comp) {
      comp_dens[n] = pi[n]*exp(normal_lpdf(y | mu[n], sigma[n]));
    }
    
    density = sum(comp_dens) + .000001;
    
    return density;
  }
  
}  

data {
  int<lower=0> N;
  vector[N] Q;
  vector<lower=0, upper=1>[N] p;
  cov_matrix[N] QCorr;
  int<lower=1> n_components;
  real m; // mean prior mean
  real<lower=0> c; // mean prior sd  
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // sample size prior sd
  // vector[n_components] alpha; // component weights prior parameter
  real a;
  real c1;
  real alph;
}

transformed data{
  real scaling_step = 1e-2;
  real rel_tol = 1e-6;
  real f_tol = 1;
  int max_steps = 500;
  vector[1] y_guess = [0.5]';
  // real alpha = 0.5;
}

parameters {
  vector[n_components] mus;
  vector<lower=0>[n_components] sigmas;
  // simplex[n_components] pi;
  real<lower=0> n;
  vector<lower=0, upper = 1>[n_components - 1] v;
}

transformed parameters {
  simplex[n_components] pi;
  real stick_remaining;

  stick_remaining = 1;
  for (k in 1:(n_components - 1)) {
    pi[k] = v[k] * stick_remaining;
    stick_remaining *= (1 - v[k]);
  }
  pi[n_components] = stick_remaining;
  
  vector<lower=0,upper=1>[N] p_hat;
  for (i in 1:N) p_hat[i] = GM_CDF(Q[i], mus, sigmas, pi, n_components);
}

model {
  v ~ beta(1, alph);
  // target += log((alpha + (a / c1))^(-1) - (alpha + ((a + 1) / c1))^(-1));
  // pi ~ dirichlet(alpha);
  mus ~ normal(m, c);
  sigmas ~ normal(0, sv);
  n ~ normal(0, nv);
  
  p_hat ~ multi_normal(p, (1/n)*QCorr);
}

generated quantities {
  real dist_samp;
  int samp_comp = categorical_rng(pi);
  
  dist_samp = normal_rng(mus[samp_comp], sigmas[samp_comp]);
  
   vector[N] p_samp;
  vector[N] Q_rep;  // Q such that GM_CDF(Q_rep[i]) ≈ p[i]
  
  p_samp = multi_normal_rng(p, (1/n)*QCorr);
  for (i in 1:N) if (p_samp[i] < 0) p_samp[i] = 1e-6;
  for (i in 1:N) if (p_samp[i] > 1) p_samp[i] = .999999;
  // p_samp[p_samp < 0] = 0;
  // p_samp[p_samp > 1] = 1;
  // for (i in 1:N) {
  //   Q_rep[i] = GM_quantile(p_samp[i], mus, sigmas, pi);
  // }
  
  for (i in 1:N) {
    Q_rep[i] = solve_powell_tol(GMInv_CDF, y_guess,
                                rel_tol, f_tol, max_steps,
                                p_samp[i], pi, mus, sigmas, N, n_components)[1];
  }
}

