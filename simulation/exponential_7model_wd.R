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

args <- commandArgs()
sample_type = args[6]
p <- as.numeric(args[7])
nind <- as.numeric(args[8])
out_s <- 5000
samples <- 5000


cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_exponential_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_exponential_ind_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_exponential_quantiles.stan')

cltnmod <- cmdstan_model(stan_file = 
                           '../stan_models/cdf_exponential_n_quantiles.stan')

ordnmod <- cmdstan_model(stan_file = 
                           '../stan_models/order_exponential_n_quantiles.stan')
mod_loc <- "../stan_models/"

#cltnmod <- paste0(mod_loc, "simple_expal_n_quantiles.stan")
#ordnmod <- paste0(mod_loc, "order_expal_n_quantiles.stan")
#cltmod <- paste0(mod_loc, "simple_expal_quantiles.stan")
#ordmod <- paste0(mod_loc, "order_expal_quantiles.stan")
#indmod <- paste0(mod_loc, "ind_simple_expal_quantiles.stan")

samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
levels <- list(
  c(.25, .5, .75),
  c(.1, .25, .5, .75, .9),
  c(.05, .1, .25, .5, .75, .9, .95),
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
ddist <- function(x) {dexp(x, lambda)}
rdist <- function(n) {rexp(n, lambda)}

models <- c("cltn", "ordn", "clt", "ord", "ind", "kern", "spline")

reps <- 1000
set.seed(92)
n <- samp_sizes[nind]
seeds <- sample(999999999, reps)
pdist <- function(x) {pexp(x, lambda)} 

distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "dplyr", "tidyr",
                                  "janitor")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
  #foreach(n = samp_sizes, .combine = rbind) %:%
  #foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
    true_params <- data.frame(variable = c("lambda", "n"),
                              truth = c(lambda, n))
    
    source("./simulation_functions.R")
    set.seed(replicate)
    samp <- rexp(n, lambda)
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs, type = 2)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = n)
    
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
    start <- Sys.time()
  fit_cltn <- stan_fit_draws(cltnmod, stan_data, sampler = sample_type)
  end <- Sys.time()
  cltn_time <- difftime(end, start, units = "min")[[1]]

  start <- Sys.time()
  fit_ordn <- stan_fit_draws(ordnmod, stan_data, sampler = sample_type)
  end <- Sys.time()
  ordn_time <- difftime(end, start, units = "min")[[1]] 
  
  start <- Sys.time()
  fit_clt <- stan_fit_draws(cltmod, stan_data, sampler = sample_type)
  end <- Sys.time()
  clt_time <- difftime(end, start, units = "min")[[1]]
  
  start <- Sys.time()
  fit_ord <- stan_fit_draws(ordmod, stan_data, sampler = sample_type)
  end <- Sys.time()
  ord_time <- difftime(end, start, units = "min")[[1]]
  
  start <- Sys.time()
  fit_ind <- stan_fit_draws(indmod,stan_data, sampler = sample_type)
  end <- Sys.time()
  ind_time <- difftime(end, start, units = "min")[[1]]

  times <- c(cltn_time, ordn_time, clt_time, ord_time, ind_time)
  
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
  
  sum_eval$time <- times
    
    
    sum_eval$model <- models[1:5]
    draws_cltn <- fit_cltn[[1]]$draws
    draws_ordn <- fit_ordn[[1]]$draws
    draws_clt <- fit_clt[[1]]$draws
    draws_ord <- fit_ord[[1]]$draws
    draws_ind <- fit_ind[[1]]$draws
    
    
    
    params <- fit_cltn[[2]]$draws(variables = "lambda", 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_lambda <- m_params[1]
    
    dcltn <- function(x) {dexp(x, cltn_lambda)}
    qcltn <- function(p) {qexp(p, cltn_lambda)}
    
    ######ordn#######
    
    params <- fit_ordn[[2]]$draws(variables = "lambda", 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    ordn_lambda <- m_params[1]
    
    dordn <- function(x) {dexp(x, ordn_lambda)}
    qordn <- function(p) {qexp(p, ordn_lambda)}
    
    
    ######clt#######
    
    params <- fit_clt[[2]]$draws(variables = "lambda", 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_lambda <- m_params[1]
    
    dclt <- function(x) {dexp(x, clt_lambda)}
    qclt <- function(p) {qexp(p, clt_lambda)}
    
    ######ord#######
    
    params <- fit_ord[[2]]$draws(variables = "lambda", 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    ord_lambda <- m_params[1]
    
    dord <- function(x) {dexp(x, ord_lambda)}
    qord <- function(p) {qexp(p, ord_lambda)}
    
    
    
    ######ind#######
    
    params <- fit_ind[[2]]$draws(variables = "lambda", 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    ind_lambda <- m_params[1]
    
    dind <- function(x) {dexp(x, ind_lambda)}
    qind <- function(p) {qexp(p, ind_lambda)}
    
    # true_draws <- rexp(length(draws_clt), mu, sigma)
    
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
    qkern <- function(p) {qkden(p, quantiles, kernel = "epanechnikov")}
    pukern <- function(x) {pdist(qkern(x))}
    rspline <- make_r_fn(probs, quantiles)
    dspline <- make_d_fn(probs, quantiles)
    rkern <- function(n) {rkden(n, quantiles, kernel = "epanechnikov")}
    dkern <- function(x) {dkden(x, quantiles, kernel = "epanechnikov")}
    pucltno <- function(x) {pdist(qcltn(x))}
    puordno <- function(x) {pdist(qordn(x))}
    puclto <- function(x) {pdist(qclt(x))}
    puordo <- function(x) {pdist(qord(x))}
    puindo <- function(x) {pdist(qind(x))}
    
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
    
    #sum_eval <- data.frame(model = mod, width_lambda = NA, width_sigma = NA,
    #			 width_n = NA, cover90_lambda = NA, cover90_sigma = NA,
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
    
    uwd1_cltno <- unit_wass_dist(pucltno, d = 1)
    uwd1_ordno <- unit_wass_dist(puordno, d = 1)
    uwd1_clto <- unit_wass_dist(puclto, d = 1)
    uwd1_ordo <- unit_wass_dist(puordo, d = 1)
    uwd1_indo <- unit_wass_dist(puindo, d = 1)
    
    
    ks_cltn <- ks.test(udraws_cltn, "punif")$statistic
    ks_ordn <- ks.test(udraws_ordn, "punif")$statistic
    ks_clt <- ks.test(udraws_clt, "punif")$statistic
    ks_ord <- ks.test(udraws_ord, "punif")$statistic
    ks_ind <- ks.test(udraws_ind, "punif")$statistic
    # ks_meta <- ks.test(udraws_meta, "punif")$statistic
    ks_spline <- ks.test(pdist(rspline(out_s)), "punif")$statistic
    ks_kern <- ks.test(rkern(out_s), "punif")$statistic
    
    
    kls <- rdist(samples)
    cltnx <- dcltn(kls)
    ordnx <- dordn(kls)
    cltx <- dclt(kls)
    ordx <- dord(kls)
    indx <- dind(kls)
    splinex <- dspline(kls)
    kernx <- dkern(kls)
    py <- ddist(kls)
    
    
    cltn_kl <- mean(log(py) - log(cltnx))
    ordn_kl <- mean(log(py) - log(ordnx))
    clt_kl <- mean(log(py) - log(cltx))
    ord_kl <- mean(log(py) - log(ordx))
    ind_kl <- mean(log(py) - log(indx))
    spline_kl <- mean(log(py) - log(splinex))
    kern_kl <- mean(log(py) - log(kernx))
    
    
    
    cltn_tv <- dens_dist(dcltn, ddist)
    ordn_tv <- dens_dist(dordn, ddist)
    clt_tv <- dens_dist(dclt, ddist)
    ord_tv <- dens_dist(dord, ddist)
    ind_tv <- dens_dist(dind, ddist)
    spline_tv <- dens_dist(dspline, ddist)
    kern_tv <- dens_dist(dkern, ddist)
    
    uwd1s <- c(uwd1_cltn, uwd1_ordn, uwd1_clt, uwd1_ord, uwd1_ind, uwd1_spline,
               uwd1_kern)
    
    uwd2s <- c(uwd2_cltn, uwd2_ordn, uwd2_clt, uwd2_ord, uwd2_ind, uwd2_spline,
               uwd2_kern)
    
    uwd1os <- c(uwd1_cltno, uwd1_ordno, uwd1_clto, uwd1_ordo, uwd1_indo, uwd1_spline,
                uwd1_kern) 
    
    kss <- c(ks_cltn, ks_ordn, ks_clt, ks_ord, ks_ind, ks_spline, ks_kern)
    
    klds <- c(cltn_kl, ordn_kl, clt_kl, ord_kl, ind_kl, spline_kl, kern_kl)
    
    tvs <- c(cltn_tv, ordn_tv, clt_tv, ord_tv, ind_tv, spline_tv, kern_tv)
    
    
    
    
    data.frame(rep = replicate, n = n, probs = p, quants = length(probs), 
               model = models, uwd1 = uwd1s, uwd1o = uwd1os, 
               uwd2 = uwd2s, ks = kss, kld = klds,
               tv = tvs) %>% 
      left_join(sum_eval, by = "model")
    
    
    
  }



write.csv(distance, paste0("sim_scores/exp_", sample_type, "/", 
                           "size", n, "_probs", length(levels[[p]]), 
                           "_scores.csv"), row.names = FALSE)



