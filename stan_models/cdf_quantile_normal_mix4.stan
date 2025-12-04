

functions{
      
  vector GMInv_CDF(vector Q, data real p, vector pi, vector mu, 
                   vector sigma, int N) {
    vector[1] z;
    real cdf;
    // this encodes the cdf from which we want to sample
    cdf = pi[1]*exp(normal_lcdf(Q | mu[1], sigma[1])) +
          + pi[2]*exp(normal_lcdf(Q | mu[2], sigma[2]))
          + pi[3]*exp(normal_lcdf(Q | mu[3], sigma[3])) 
          + pi[4]*exp(normal_lcdf(Q | mu[4], sigma[4]))
          // + pi[5]*exp(normal_lcdf(Q | mu[5], sigma[5]))
          ;
    z[1] = cdf - p;
    return z;
  }
  
  real GM_CDF(real y, vector mu, vector sigma, vector pi) {
    real prob;
    
    prob = pi[1]*normal_cdf(y | mu[1], sigma[1])
           + pi[2]*normal_cdf(y | mu[2], sigma[2])
           + pi[3]*normal_cdf(y | mu[3], sigma[3])
           + pi[4]*normal_cdf(y | mu[4], sigma[4])
           // + pi[5]*normal_cdf(y | mu[5], sigma[5])
           ;
           
    return prob;
  }
  
  real GM_PDF(vector mu, vector sigma, vector pi, real y) {
    
    real density;
    
    density = pi[1]*exp(normal_lpdf(y | mu[1], sigma[1])) +
              + pi[2]*exp(normal_lpdf(y | mu[2], sigma[2]))
              + pi[3]*exp(normal_lpdf(y | mu[3], sigma[3]))
              + pi[4]*exp(normal_lpdf(y | mu[4], sigma[4]))
              // + pi[5]*exp(normal_lpdf(y | mu[5], sigma[5]))
              + .000001
              ;
    
    return density;
  }
  
  
  real GM_quantile(real p, vector mu, vector sigma, vector pi) {
    // Quantile by bisection
    real left = -10;
    real right =  10;
    real mid;
    real cdf_mid;

    for (n in 1:60) {                     // 2^-60 ≈ 1e-18 precision
      mid = 0.5 * (left + right);
      cdf_mid = GM_CDF(mid, mu, sigma, pi);
      if (cdf_mid > p)
        right = mid;
      else
        left = mid;
    }
    return mid;
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
  vector[n_components] alpha; // component weights prior parameter
}

transformed data{
  real scaling_step = 1e-2;
  real rel_tol = 1e-6;
  real f_tol = 1;
  int max_steps = 500;
  vector[1] y_guess = [0.5]';
}

parameters {
  vector[n_components] mus;
  vector<lower=0>[n_components] sigmas;
  simplex[n_components] pi;
  real<lower=0> n;
}

transformed parameters {
  vector<lower=0,upper=1>[N] p_hat;
  for (i in 1:N) p_hat[i] = GM_CDF(Q[i], mus, sigmas, pi);
}

model {
  pi ~ dirichlet(alpha);
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
                                p_samp[i], pi, mus, sigmas, N)[1];
  }
  
}

