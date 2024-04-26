source("./simulation_functions.R")
library(cmdstanr)
#library(distfromq)
#library(evmix)
library(parallel)
library(doParallel)
library(doMC)
library(distfromq)
library(evmix)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

cltmod <- cmdstan_model(stan_file = 
                       '../stan_models/normal_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/ind_normal_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_quantiles.stan')

cltnmod <- cmdstan_model(stan_file = 
                           '../stan_models/normal_n_quantiles.stan')

ordnmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_n_quantiles.stan')
mod_loc <- "../stan_models/"

#cltnmod <- paste0(mod_loc, "simple_normal_n_quantiles.stan")
#ordnmod <- paste0(mod_loc, "order_normal_n_quantiles.stan")
#cltmod <- paste0(mod_loc, "simple_normal_quantiles.stan")
#ordmod <- paste0(mod_loc, "order_normal_quantiles.stan")
#indmod <- paste0(mod_loc, "ind_simple_normal_quantiles.stan")

samp_sizes <- c(50, 100, 500, 1000, 5000)
levels <- list(
	       c(.25, .5, .75),
	       c(.25, .75),
	       c(.025, .25, .75, .975),
	       c(.025, .25, .5, .75, .975),

	       c(.025, seq(.05, .95, by = .05), .975),
	       c(.01, .025, seq(.05, .95, by = .05), .975, .99))

mu <- 1
sigma <- 2.2
qtrue <- function(p) {qnorm(p, mu, sigma)}

models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern")

# probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
probs <- seq(0.05, .95, by = 0.05)
# probs <- c(.025, .25, .5, .75, .975)
text <- qnorm(.99, mu, sigma)
reps <- 10

distance <- foreach(replicate = 1:reps,
		    .packages = c("cmdstanr", "evmix", "distfromq")
		    ,.errorhandling = "remove"
		    ,.combine = rbind) %:%
		foreach(n = samp_sizes, .combine = rbind) %:%
			foreach(p = 1:length(levels), .combine = rbind) %dopar% {
 
#replicate = 2
#n = 100
#p = 3
source("./simulation_functions.R")
  samp <- rnorm(n, mu, sigma)
  probs <- levels[[p]]
  quantiles <- quantile(samp, probs)
  
  data <- data.frame(quantile = quantiles, prob = probs)
  stan_data <- make_stan_data(data, size = n)
  
  #fit models
  draws_cltn <- stan_fit_draws(cltnmod, stan_data)
  draws_ordn <- stan_fit_draws(ordnmod, stan_data)
  draws_clt <- stan_fit_draws(cltmod, stan_data)
  draws_ord <- stan_fit_draws(ordmod, stan_data)
  draws_ind <- stan_fit_draws(indmod,stan_data)
  
   
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
  wd1_ordn <- wass_dist(qordn, qtrue)
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
  
  
  #print("rep, n, p, length(probs), models, wd1s, wd2s") 
  
  
  
  
  #unit draws
  udraws_cltn <- pnorm(draws_cltn, mu, sigma)
  udraws_ordn <- pnorm(draws_ordn, mu, sigma)
  udraws_clt <- pnorm(draws_clt, mu, sigma)
  udraws_ord <- pnorm(draws_ord, mu, sigma)
  udraws_ind <- pnorm(draws_ind, mu, sigma)
  
  
  #make unit ecdfs
  pucltn <- function(x) {ecdf(udraws_cltn)(x)}
  puordn <- function(x) {ecdf(udraws_ordn)(x)}
  puclt <- function(x) {ecdf(udraws_clt)(x)}
  puord <- function(x) {ecdf(udraws_ord)(x)}
  puind <- function(x) {ecdf(udraws_ind)(x)}
  qspline <- make_q_fn(probs, quantiles)
  puspline <- function(x) {pnorm(qspline(x), mu, sigma)}
  qkern <- function(p) {qkden(p, quantiles, kernel = "triangular")}
  pukern <- function(x) {pnorm(qkern(x), mu, sigma)}
  
  
  uwd1_cltn <- unit_wass_dist(pucltn, d = 1)
  uwd1_ordn <- unit_wass_dist(puordn, d = 1)
  uwd1_clt <- unit_wass_dist(puclt, d = 1)
  uwd1_ord <- unit_wass_dist(puord, d = 1)
  uwd1_ind <- unit_wass_dist(puind, d = 1)
  uwd1_spline <- unit_wass_dist(puspline, d = 1)
  uwd1_kern <- unit_wass_dist(pukern, d = 1)
  
  uwd1s <- c(uwd1_cltn, uwd1_ordn, uwd1_clt, uwd1_ord, uwd1_ind, uwd1_spline,
             uwd1_kern)
  
  data.frame(rep = replicate, n = n, probs = p, quants = length(probs), model = models, 
             wd1 = wd1s, wd2 = wd2s, uwd1 = uwd1s)
  
  

}


write.csv(distance, "test.csv", row.names = FALSE)



