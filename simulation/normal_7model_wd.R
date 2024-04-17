library(cmdstanr)
library(distfromq)
library(evmix)

cltmod <- cmdstan_model(stan_file = 
                       './stan_models/simple_normal_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          './stan_models/ind_simple_normal_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          './stan_models/order_normal_quantiles.stan')

cltnmod <- cmdstan_model(stan_file = 
                           './stan_models/simple_normal_n_quantiles.stan')

ordnmod <- cmdstan_model(stan_file = 
                          './stan_models/order_normal_n_quantiles.stan')

n <- 100
mu <- 1
sigma <- 2.2
qtrue <- function(p) {qnorm(p, mu, sigma)}

models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern")

# probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
probs <- seq(0.05, .95, by = 0.05)
# probs <- c(.025, .25, .5, .75, .975)
text <- qnorm(.99, mu, sigma)
reps <- 1000
rep <- 1
norm_better <- c()
ordnorm_better <- c()
cltnorm_dist <- c()
indnorm_dist <- c()
ordnorm_dist <- c()
cltnnorm_dist <- c()
ordnnorm_dist <- c()
repeat {

  samp <- rnorm(n, mu, sigma)
  quantiles <- quantile(samp, probs)
  
  data <- data.frame(quantile = quantiles, prob = probs)
  stan_data <- make_stan_data(data)
  
  
  #fit models
  draws_cltn <- stan_fit_draws(cltnmod, stan_data)
  draws_ordn <- stan_fit_draws(ordnmod, stan_data)
  draws_clt <- stan_fit_draws(cltmod, stan_data)
  draws_ord <- stan_fit_draws(ordmod, stan_data)
  draws_ind <- stan_fit_draws(indmod, stan_data)
  
  
  #make quantile functions
  qcltn <- function(p) {quantile(draws_cltn, p)}
  qordn <- function(p) {quantile(draws_ordn, p)}
  qclt <- function(p) {quantile(draws_clt, p)}
  qord <- function(p) {quantile(draws_ord, p)}
  qind <- function(p) {quantile(draws_ind, p)}
  #spline and kernel models go directly to quantile function
  qspline <- make_q_fn(probs, quantiles)
  qkern <- function(p) {qkden(p, quantiles, kernel = "triangular")}
  
  
  #calculate 1wd
  wd1_cltn <- wass_dist(qcltn, qtrue)
  wd1_ordn <- wass_dist(ordn, qtrue)
  wd1_clt <- wass_dist(qclt, qtrue)
  wd1_ord <- wass_dist(qord, qtrue)
  wd1_ind <- wass_dist(qind, qtrue)
  wd1_spline <- wass_dist(qspline, qtrue)
  wd1_kern <- wass_dist(qkern, qtrue)
  
  
  #calculate 2wd
  wd2_cltn <- wass_dist(qcltn, qtrue, d = 2)
  wd2_ordn <- wass_dist(qordn, qtrue, d = 2)
  wd2_clt <- wass_dist(qclt, qtrue, d = 2)
  wd2_ord <- wass_dist(qord, qtrue, d = 2)
  wd2_ind <- wass_dist(qind, qtrue, d = 2)
  wd2_spline <- wass_dist(qspline, qtrue, d = 2)
  wd2_kern <- wass_dist(qkern, qtrue, d = 2)
  
  wd1s <- c(wd1_cltn, wd1_ordn, wd1_clt, wd1_ord, wd1_ind, wd1_spline, wd1_kern)
  wd2s <- c(wd2_cltn, wd2_ordn, wd2_clt, wd2_ord, wd2_ind, wd2_spline, wd2_kern)
  
  
  
  data.frame(rep = m, n = n, quants = length(probs), model = models, 
             wd1 = wd1s, wd2 = wd2s)
  
  rep <- rep + 1
  if (rep > reps) {break}

}






