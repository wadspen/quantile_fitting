


functions {
  
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
  
  
  real GM_quantile(real p, vector mu, vector sigma, vector pi, int num_comp) {
    // Quantile by bisection
    real left = -10;
    real right =  10;
    real mid;
    real cdf_mid;

    for (n in 1:60) {                     // 2^-60 ≈ 1e-18 precision
      mid = 0.5 * (left + right);
      cdf_mid = GM_CDF(mid, mu, sigma, pi, num_comp);
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
  int<lower=1> n_components;
  real m;
  real<lower=0> c; 
  real<lower=0> sv; // sd prior sd
  real<lower=0> nv; // n prior sd
  vector[n_components] alpha;
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
  real<lower=0> n;
  // vector<lower=0,upper=1>[n_components] pi;
  simplex[n_components] pi;
}

transformed parameters {
  vector[N] U;
  // vector[n_components] pit;
  
  // pit = softmax(pi);
  
  for (i in 1:N)
  
    U[i] = GM_CDF(Q[i], mus, sigmas, pi, n_components);
    // U[i] = pi[1]*normal_cdf(Q[i]| mus[1], sigmas[1]) + 
    //        pi[2]*normal_cdf(Q[i]| mus[2], sigmas[2]) +
    //        pi[3]*normal_cdf(Q[i]| mus[3], sigmas[3])
    //        + pi[4]*normal_cdf(Q[i]| mus[4], sigmas[4])
    //        //+ pi[5]*normal_cdf(Q[i] | mus[5], sigmas[5])
    //        ;
}

model {
  // pi ~ normal(4,5);
  pi ~ dirichlet(alpha);
  mus ~ normal(m, c);
  sigmas ~ normal(0, sv);
  n ~ normal(0, nv);

  
  target += orderstatistics(n, N, p, U);
  
  
  for (i in 1:N) {
    vector[n_components] lps;
    for (k in 1:n_components) {
      lps[k] = log(pi[k]) + normal_lpdf(Q[i] | mus[k], sigmas[k]);
    }
    target += log_sum_exp(lps);
    
  }
  
  // for (i in 1:N)
  //   target += log_sum_exp({log(pi[1]) + normal_lpdf(Q[i] | mus[1], sigmas[1]),
  //                         log(pi[2]) + normal_lpdf(Q[i] | mus[2], sigmas[2]),
  //                         log(pi[3]) + normal_lpdf(Q[i] | mus[3], sigmas[3]),
  //                         log(pi[4]) + normal_lpdf(Q[i] | mus[4], sigmas[4])
  //                         //,log(pi[5]) + normal_lpdf(Q[i] | mus[5], sigmas[5])
  //                         });
    
}

generated quantities {
  vector[N] pred_q;
  real dist_samp;
  int samp_comp = categorical_rng(pi);
  dist_samp = normal_rng(mus[samp_comp], sigmas[samp_comp]);
  
  vector<lower=0, upper=1>[N] p_samp;
  vector[N] Q_rep;  // Q such that GM_CDF(Q_rep[i]) ≈ p[i]
  
  for (i in 1:N) p_samp[i] = beta_rng(n*p[i], n - p[i]*n + 1);
  
  // for (i in 1:N) {
  //   Q_rep[i] = GM_quantile(p_samp[i], mus, sigmas, pi);
  // }
  
  for (i in 1:N) {
    Q_rep[i] = solve_powell_tol(GMInv_CDF, y_guess,
                                rel_tol, f_tol, max_steps,
                                p_samp[i], pi, mus, sigmas, N, n_components)[1];
  }
}

