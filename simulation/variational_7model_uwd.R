source("./simulation_functions.R")
library(cmdstanr)
library(distr)
library(distfromq)
library(evmix)
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
dist <- args[6]

cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/normal_mix4_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/normal_mix4_ind_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_mix4_quantiles.stan')

cltnmod <- cmdstan_model(stan_file =
                           '../stan_models/normal_n_mix4_quantiles.stan')

ordnmod <- cmdstan_model(stan_file =
                           '../stan_models/order_normal_n_mix4_quantiles.stan')

mod_loc <- "../stan_models/"

#cltnmod <- paste0(mod_loc, "simple_normal_n_quantiles.stan")
#ordnmod <- paste0(mod_loc, "order_normal_n_quantiles.stan")
#cltmod <- paste0(mod_loc, "simple_normal_quantiles.stan")
#ordmod <- paste0(mod_loc, "order_normal_quantiles.stan")
#indmod <- paste0(mod_loc, "ind_simple_normal_quantiles.stan")

samp_sizes <- c(25, 50, 100, 500, 1000, 2000)
levels <- list(
  c(.4, .5, .6),
  c(.05, .4, .5, .6, .95),
  c(.3, .4, .5, .6, .7),
  c(.2, .3, .4, .5, .6, .7, .8),
  c(.01, .1, .2, .25, .5, .75, .8, .9, .99),
  c(.1, .2, .3, .4, .5, .6, .7, .8, .9),
  c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99),
  seq(.15, .85, by = .05),
  seq(.1, .9, by = .05),
  seq(.05, .95, by = .05),
  c(.025, seq(.05, .95, by = .05), .975),
  c(.01, .025, seq(.05, .95, by = .05), .975, .99))

tails <- c(0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1)
xi <- 1
omega <- 2.2
alpha <- 12
qtrue <- function(p) {qsn(p, xi, omega, alpha)}

models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern")

# probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
probs <- seq(0.05, .95, by = 0.05)
reps <- 1500

distance <- foreach(rep = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                  "VGAM", "distr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %:%
  foreach(n = samp_sizes, .combine = rbind) %:%
  foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
    #replicate = 2
    #n = 100
    #p = 3
    source("./simulation_functions.R")
    # samp <- rsn(n, xi, omega, alpha)
    # samp <- r(mixdist)(n)
    if (dist == "norm") {
      samp <- rnorm(n)
      pdist <- function(x) {pnorm(x)}
    } else if (dist == "evd") {
      samp <- revd(n)
      pdist <- function(x) {pevd(x)}
    } else if (dist == "lp") {
      samp <- rlaplace(n)
      pdist <- function(x) {plaplace(x)}
    } else if (dist == "gmix") {
      pars <- data.frame(mu = c(-.7, 1.2),
                         sigma = c(.9, .6),
                         weight = c(.35, .65))
      mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
      samp <- r(mdist)(n)
      pdist <- function(x) {p(mdist)(x)}
    }
    
    tail <- tails[p]
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = n, comps = 3)
    
    #fit models
    draws_cltn <- stan_fit_draws(cltnmod, stan_data, 
                                 sampler = "variational",
                                 elbo = 400, grad = 25,
                                 out_s = 2000)
    draws_ordn <- stan_fit_draws(ordnmod, stan_data, 
                                 sampler = "variational",
                                 elbo = 400, grad = 25,
                                 out_s = 2000)
    draws_clt <- stan_fit_draws(cltmod, stan_data,
                                sampler = "variational",
                                elbo = 400, grad = 25,
                                out_s = 2000)
    draws_ord <- stan_fit_draws(ordmod, stan_data,
                                sampler = "variational",
                                elbo = 400, grad = 25,
                                out_s = 2000)
    draws_ind <- stan_fit_draws(indmod,stan_data,
                                sampler = "variational",
                                elbo = 400, grad = 25,
                                out_s = 2000)
    
    
    #unit draws
    udraws_cltn <- pdist(draws_cltn)
    udraws_ordn <- pdist(draws_ordn)
    udraws_clt <- pdist(draws_clt)
    udraws_ord <- pdist(draws_ord)
    udraws_ind <- pdist(draws_ind)
    
    
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
    
    data.frame(rep = rep, n = n, probs = p, quants = length(probs), 
               model = models, 
               uwd1 = uwd1s, uwd2 = uwd2s)
    
    
    
  }


write.csv(distance, paste0(dist, "_test.csv"), row.names = FALSE)



