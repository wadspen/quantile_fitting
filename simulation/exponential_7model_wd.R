source("./simulation_functions.R")
library(cmdstanr)
library(distfromq)
library(evmix)
library(bayesplot)
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/exponential_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/exponential_ind_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_exponential_quantiles.stan')

cltnmod <- cmdstan_model(stan_file = 
                           '../stan_models/exponential_n_quantiles.stan')

ordnmod <- cmdstan_model(stan_file = 
                           '../stan_models/order_exponential_n_quantiles.stan')
mod_loc <- "../stan_models/"

#cltnmod <- paste0(mod_loc, "simple_normal_n_quantiles.stan")
#ordnmod <- paste0(mod_loc, "order_normal_n_quantiles.stan")
#cltmod <- paste0(mod_loc, "simple_normal_quantiles.stan")
#ordmod <- paste0(mod_loc, "order_normal_quantiles.stan")
#indmod <- paste0(mod_loc, "ind_simple_normal_quantiles.stan")

samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
levels <- list(
  seq(.1, .9, by = .1),
  c(.05, seq(.1, .9, by = .1), .95),
  c(.025, .05, seq(.1, .9, by = .1), .95, .975),
  c(.01, .025, .05, seq(.1, .9, by = .1), .95, .975, .99),
  seq(.05, .95, by = .05),
  c(.025, seq(.05, .95, by = .05), .975),
  c(.01, .025, seq(.05, .95, by = .05), .975, .99),
  seq(.01, .99, by = .02)	       
)



lambda <- .4
qtrue <- function(p) {qexp(p, lambda)}

models <- c("cltn", "ordn", "clt", "ord", "ind", "kern", "spline")

reps <- 1000
set.seed(92)

seeds <- sample(999999999, reps)
pdist <- function(x) {pexp(x, lambda)} 

distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "dplyr", "tidyr",
                                  "janitor")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %:%
  foreach(n = samp_sizes, .combine = rbind) %:%
  foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
    true_params <- data.frame(variable = c("lambda", "n"),
                              truth = c(lambda, n))
    
    source("./simulation_functions.R")
    set.seed(replicate)
    samp <- rexp(n, lambda)
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs, type = 2)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = n, cor = "exp")
    
    #fit models
    #if (mod %in% c("cltn", "ordn", "clt", "ord", "ind")) {
    #if (mod == "cltn") {
    #	fit <- stan_fit_draws(cltnmod, stan_data)
    #} else if (mod == "ordn") {
    #	fit <- stan_fit_draws(ordnmod, stan_data)
    #} else if (mod == "clt") {
    #	fit <- stan_fit_draws(cltmod, stan_data)
    #} else if (mod == "ord") {
    #	fit <- stan_fit_draws(ordmod, stan_data)
    #} else if (mod == "ind") {
    #	fit <- stan_fit_draws(indmod, stan_data)
    #}
    
    # samp <- cltmod$sample(stan_data, init = list(list(lambda = .399),
    #                                               list(lambda = .399),
    #                                               list(lambda = .399),
    #                                               list(lambda = .399)))
    # samp$summary()
    # mcmc_trace(samp$draws(), pars = "lambda")
    #sum <- fit$summary
    #sum_eval <- eval_sum(sum, true_params)	
    fit_cltn <- stan_fit_draws(cltnmod, stan_data)
    fit_ordn <- stan_fit_draws(ordnmod, stan_data)
    fit_clt <- stan_fit_draws(cltmod, stan_data)
    fit_ord <- stan_fit_draws(ordmod, stan_data)
    fit_ind <- stan_fit_draws(indmod,stan_data)
    
    sum_cltn <- fit_cltn[[2]]$summary()
    sum_ordn <- fit_ordn[[2]]$summary()
    sum_clt <- fit_clt[[2]]$summary()
    sum_ord <- fit_ord[[2]]$summary()
    sum_ind <- fit_ind[[2]]$summary()
    
    sum_eval <- rbind(eval_sum(sum_cltn, true_params),
                      eval_sum(sum_ordn, true_params),
                      eval_sum(sum_clt, true_params),
                      eval_sum(sum_ord, true_params),
                      eval_sum(sum_ind, true_params)
    )
    
    
    sum_eval$model <- models[1:5]
    draws_cltn <- fit_cltn[[1]]$draws
    draws_ordn <- fit_ordn[[1]]$draws
    draws_clt <- fit_clt[[1]]$draws
    draws_ord <- fit_ord[[1]]$draws
    draws_ind <- fit_ind[[1]]$draws
    
    # true_draws <- rnorm(length(draws_clt), mu, sigma)
    
    # kldiv(draws_ind, true_draws)
    #unit draws
    udraws_cltn <- pdist(draws_cltn)
    udraws_ordn <- pdist(draws_ordn)
    udraws_clt <- pdist(draws_clt)
    udraws_ord <- pdist(draws_ord)
    udraws_ind <- pdist(draws_ind)
    
    #udraws <- pdist(draws)
    
    #make unit ecdfs
    pucltn <- function(x) {ecdf(udraws_cltn)(x)}
    puordn <- function(x) {ecdf(udraws_ordn)(x)}
    puclt <- function(x) {ecdf(udraws_clt)(x)}
    puord <- function(x) {ecdf(udraws_ord)(x)}
    puind <- function(x) {ecdf(udraws_ind)(x)}
    qspline <- make_q_fn(probs, quantiles)
    puspline <- function(x) {pdist(qspline(x))}
    qkern <- function(p) {qkden(p, quantiles, kernel = "gaussian")}
    pukern <- function(x) {pdist(qkern(x))}
    
    #pu <- function(x) {ecdf(udraws)(x)}
    #uwd1 <- unit_wass_dist(pu, d = 1)
    #uwd2 <- unit_wass_dist(pu, d = 2)
    #} else if (mod %in% c("kern", "spline")) {
    #	  if (mod == "kern") {
    #		  qu <- function(p) {qkden(p, quantiles, kernel = "gaussian")}
    #        } else if (mod == "spline") {
    #		  qu <- make_q_fn(probs, quantiles)
    #	  }
    
    #pu <- function(x) {pdist(qu(x))}
    #uwd1 <- unit_wass_dist(pu, d = 1)
    #uwd2 <- unit_wass_dist(pu, d = 2)
    
    #sum_eval <- data.frame(model = mod, width_mu = NA, width_sigma = NA,
    #			 width_n = NA, cover90_mu = NA, cover90_sigma = NA,
    #			 cover90_n = NA)
    #}
    
    # draws <- draws_ord
    # true_draws <- rexp(length(draws_ord), lambda)
    # 
    # kldiv(draws, true_draws)
    
    uwd1_cltn <- unit_wass_dist(pucltn, d = 1)
    uwd1_ordn <- unit_wass_dist(puordn, d = 1)
    uwd1_clt <- unit_wass_dist(puclt, d = 1)
    uwd1_ord <- unit_wass_dist(puord, d = 1)
    uwd1_ind <- unit_wass_dist(puind, d = 1)
    uwd1_spline <- unit_wass_dist(puspline, d = 1)
    uwd1_kern <- unit_wass_dist(pukern, d = 1)
    
    uwd2_cltn <- unit_wass_dist(pucltn, d = 2)
    uwd2_ordn <- unit_wass_dist(puordn, d = 2)
    uwd2_clt <- unit_wass_dist(puclt, d = 2)
    uwd2_ord <- unit_wass_dist(puord, d = 2)
    uwd2_ind <- unit_wass_dist(puind, d = 2)
    uwd2_spline <- unit_wass_dist(puspline, d = 2)
    uwd2_kern <- unit_wass_dist(pukern, d = 2)
    
    uwd1s <- c(uwd1_cltn, uwd1_ordn, uwd1_clt, uwd1_ord, uwd1_ind, uwd1_spline,
               uwd1_kern)
    
    uwd2s <- c(uwd2_cltn, uwd2_ordn, uwd2_clt, uwd2_ord, uwd2_ind, uwd2_spline,
               uwd2_kern)
    
    data.frame(rep = replicate, n = n, probs = p, quants = length(probs), 
               model = models, uwd1 = uwd1s, uwd2 = uwd2s) %>% 
      left_join(sum_eval, by = "model")
    
    
    
  }


write.csv(distance, "results2.csv", row.names = FALSE)



