source("./simulation_functions.R")
setwd("~/quantile_fitting/simulation/")
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
                          '../stan_models/normal_t_mix4_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/normal_t_mix4_ind_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_mix4_quantiles.stan')

cltnmod <- cmdstan_model(stan_file =
                           '../stan_models/normal_n_mix4_quantiles.stan')

ordnmod <- cmdstan_model(stan_file =
                           '../stan_models/order_normal_n_mix4_quantiles.stan')

metamod <- cmdstan_model(stan_file = 
                           "../stan_models/metanorm_quantiles.stan")

normmod <- cmdstan_model(stan_file =
			   "../stan_models/normal_quantiles.stan")

mod_loc <- "../stan_models/"

#cltnmod <- paste0(mod_loc, "simple_normal_n_quantiles.stan")
#ordnmod <- paste0(mod_loc, "order_normal_n_quantiles.stan")
#cltmod <- paste0(mod_loc, "simple_normal_quantiles.stan")
#ordmod <- paste0(mod_loc, "order_normal_quantiles.stan")
#indmod <- paste0(mod_loc, "ind_simple_normal_quantiles.stan")

samp_sizes <- c(50, 150, 500, 1000)
levels <- list(
  #c(.2, .3, .4),
  c(.25, .5, .75),
  #c(.05, .4, .5, .6, .95),
  c(.025, .25, .5, .75, .975),
  #c(.2, .3, .4, .5, .6, .7, .8),
  #c(.01, .1, .2, .25, .5, .75, .8, .9, .99),
  c(.1, .2, .3, .4, .5, .6, .7, .8, .9),
  #c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99),
  #seq(.15, .85, by = .05),
  #seq(.1, .9, by = .05),
  #seq(.05, .95, by = .05),
  c(.025, seq(.05, .95, by = .05), .975))#,
#c(.01, .025, seq(.05, .95, by = .05), .975, .99))

#tails <- c(.5, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1)
tails <- c(0, 0, 1)


tails <- c(0, 1, 0, 1)
xi <- 1
omega <- 2.2
alpha <- 12
qtrue <- function(p) {qsn(p, xi, omega, alpha)}

models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern", "meta", "norm")

# probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
probs <- seq(0.05, .95, by = 0.05)
reps <- 1000

distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                  "VGAM", "distr", "dplyr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %:%
  foreach(n = samp_sizes, .combine = rbind) %:%
  foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
    
    source("./simulation_functions.R")
    
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
      pars <- data.frame(mu = c(-1, 1.2),
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
    stan_data <- make_stan_data(data, size = n, comps = 5)
    
    
    #fit models
    fit_cltn <- stan_fit_draws(cltnmod, stan_data, 
                                 sampler = "variational",
                                 elbo = 300, grad = 10,
                                 out_s = 2000, refresh = 100)
    fit_ordn <- stan_fit_draws(ordnmod, stan_data, 
                                 sampler = "variational",
                                 elbo = 300, grad = 10,
                                 out_s = 2000)
    fit_clt <- stan_fit_draws(cltmod, stan_data,
                               sampler = "MCMC", refresh = 10, burn = 2000,
                              samp = 2000,
                               elbo = 300, grad = 10,
                               out_s = 10000)
    fit_ord <- stan_fit_draws(ordmod, stan_data,
                                sampler = "MCMC", refresh = 10,
                                burn = 2000, samp = 2000,
                                elbo = 300, grad = 10,
                                out_s = 10000)
    
    fit_ind <- stan_fit_draws(indmod,stan_data,
                                sampler = "MCMC", burn = 2000, samp = 2000,
                                refresh = 10,
                                elbo = 300, grad = 10,
                                out_s = 2000)
    
    fit_meta <- stan_fit_draws(metamod,stan_data,
                                sampler = "variational",
                                elbo = 300, grad = 10,
                                out_s = 2000)

    fit_norm <- stan_fit_draws(normmod,stan_data,
			       sampler = "variational",
			       elbo = 300, grad = 10,
			       out_s = 2000)
    
    
    #unit draws
    udraws_cltn <- pdist(fit_cltn$draws)
    udraws_ordn <- pdist(fit_ordn$draws)
    udraws_clt <- pdist(fit_clt$draws)
    udraws_ord <- pdist(fit_ord$draws)
    udraws_ind <- pdist(fit_ind$draws)
    udraws_meta <- pdist(fit_meta$draws)
    udraws_norm <- pdist(fit_norm$draws)
    
    
    #make unit ecdfs
    pucltn <- function(x) {ecdf(udraws_cltn)(x)}
    puordn <- function(x) {ecdf(udraws_ordn)(x)}
    puclt <- function(x) {ecdf(udraws_clt)(x)}
    puord <- function(x) {ecdf(udraws_ord)(x)}
    puind <- function(x) {ecdf(udraws_ind)(x)}
    pumeta <- function(x) {ecdf(udraws_meta)(x)}
    punorm <- function(x) {ecdf(udraws_norm)(x)}
    qspline <- make_q_fn(probs, quantiles)
    puspline <- function(x) {pdist(qspline(x))}
    qkern <- function(p) {qkden(p, quantiles, kernel = "epanechnikov")}
    pukern <- function(x) {pdist(qkern(x))}
    
    
    uwd1_cltn <- unit_wass_dist(pucltn, d = 1)
    uwd1_ordn <- unit_wass_dist(puordn, d = 1)
    uwd1_clt <- unit_wass_dist(puclt, d = 1)
    uwd1_ord <- unit_wass_dist(puord, d = 1)
    uwd1_ind <- unit_wass_dist(puind, d = 1)
    uwd1_meta <- unit_wass_dist(pumeta, d = 1)
    uwd1_norm <- unit_wass_dist(punorm, d = 1)
    uwd1_spline <- unit_wass_dist(puspline, d = 1)
    uwd1_kern <- unit_wass_dist(pukern, d = 1)
    
    uwd2_cltn <- unit_wass_dist(pucltn, d = 2)
    uwd2_ordn <- unit_wass_dist(puordn, d = 2)
    uwd2_clt <- unit_wass_dist(puclt, d = 2)
    uwd2_ord <- unit_wass_dist(puord, d = 2)
    uwd2_ind <- unit_wass_dist(puind, d = 2)
    uwd2_meta <- unit_wass_dist(pumeta, d = 2)
    uwd2_norm <- unit_wass_dist(punorm, d = 2)
    uwd2_spline <- unit_wass_dist(puspline, d = 2)
    uwd2_kern <- unit_wass_dist(pukern, d = 2)
    
    uwd1s <- c(uwd1_cltn, uwd1_ordn, uwd1_clt, uwd1_ord, uwd1_ind, uwd1_spline,
               uwd1_kern, uwd1_meta, uwd1_norm)
    uwd2s <- c(uwd2_cltn, uwd2_ordn, uwd2_clt, uwd2_ord, uwd2_ind, uwd2_spline,
               uwd2_kern, uwd2_meta, uwd1_norm)
    
     
    data.frame(rep = replicate, n = n, probs = p, quants = length(probs), 
               tail = tails[p], model = models, uwd1 = uwd1s, uwd2 = uwd2s)

    
    
    
    
  }


write.csv(distance, paste0(dist, "_test.csv"), row.names = FALSE)



