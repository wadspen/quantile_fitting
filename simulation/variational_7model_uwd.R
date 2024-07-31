source("./simulation_functions.R")
#setwd("./simulation/")
library(cmdstanr)
library(distr)
library(distfromq)
library(evmix)
library(dplyr)
library(VGAM)
library(EnvStats)
library(tidyr)
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
p <- as.numeric(args[7])
nind <- as.numeric(args[8])

print(dist); print(p); print(nind)


cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_quantile_normal_mix5.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_ind_quantile_normal_mix5.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_mix5_quantiles.stan')

cltnmod <- cmdstan_model(stan_file =
                           '../stan_models/cdf_quantile_normal_n_mix5.stan')

ordnmod <- cmdstan_model(stan_file =
                           '../stan_models/order_normal_n_mix5_quantiles.stan')

metamod <- cmdstan_model(stan_file = 
                           "../stan_models/metanorm_quantiles.stan")

#normmod <- cmdstan_model(stan_file =
#			   "../stan_models/normal_quantiles.stan")

mod_loc <- "../stan_models/"

burn <- 6000
samples <- 7000
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

models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern", "meta")


reps <- 50
n <- samp_sizes[nind]


#p <- 2
#n <- 500
#dist <- "norm"
#replicate <- 2
distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                  "VGAM", "distr", "dplyr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
  #foreach(n = samp_sizes, .combine = rbind) %:%
  #foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
   
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
    
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs, type = 2)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = n, comps = 4)
    
    
    #fit models
    fit_cltn <- stan_fit_draws(cltnmod, stan_data, 
                                sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
	saveRDS(list(draws = fit_cltn[[2]],samp = samp), paste0("sim_draws5/", dist, "_cltn_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
   print("gets here") 
    fit_ordn <- stan_fit_draws(ordnmod, stan_data, 
                                sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 

   	saveRDS(list(draws = fit_ordn[[2]],samp = samp), paste0("sim_draws5/", dist, "_ordn_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
    fit_clt <- stan_fit_draws(cltmod, stan_data,
                               sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 

   	saveRDS(list(draws = fit_clt[[2]],samp = samp), paste0("sim_draws5/", dist, "_clt_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
    fit_ord <- stan_fit_draws(ordmod, stan_data,
                                sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
    
   	saveRDS(list(draws = fit_ord[[2]],samp = samp), paste0("sim_draws5/", dist, "_ord_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
    fit_ind <- stan_fit_draws(indmod,stan_data,
                                sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000)
   
	saveRDS(list(draws = fit_ind[[2]],samp = samp), paste0("sim_draws5/", dist, "_ind_rep", replicate, "_size", n, "_probs", length(probs), ".rds"))
    fit_meta <- stan_fit_draws(metamod,stan_data,
                               sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000)
	saveRDS(list(draws = fit_meta[[2]],samp = samp), paste0("sim_draws5/", dist, "_meta_rep", replicate, "_size", n, "_probs", length(probs), ".rds"))
#    fit_norm <- stan_fit_draws(normmod,stan_data,
#			       sampler = "variational",
#			       elbo = 300, grad = 10,
#			       out_s = 2000)
    
    


    fit_cltn <- fit_cltn[[1]]
    fit_ordn <- fit_ordn[[1]]
    fit_clt <- fit_clt[[1]]
    fit_ord <- fit_ord[[1]]
    fit_ind <- fit_ind[[1]]
    fit_meta <- fit_meta[[1]]



    #unit draws
    udraws_cltn <- pdist(fit_cltn$draws)
    udraws_ordn <- pdist(fit_ordn$draws)
    udraws_clt <- pdist(fit_clt$draws)
    udraws_ord <- pdist(fit_ord$draws)
    udraws_ind <- pdist(fit_ind$draws)
    udraws_meta <- pdist(fit_meta$draws)
    #udraws_norm <- pdist(fit_norm$draws)
    
    
    #make unit ecdfs
    pucltn <- function(x) {ecdf(udraws_cltn)(x)}
    puordn <- function(x) {ecdf(udraws_ordn)(x)}
    puclt <- function(x) {ecdf(udraws_clt)(x)}
    puord <- function(x) {ecdf(udraws_ord)(x)}
    puind <- function(x) {ecdf(udraws_ind)(x)}
    pumeta <- function(x) {ecdf(udraws_meta)(x)}
    #punorm <- function(x) {ecdf(udraws_norm)(x)}
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
    #uwd1_norm <- unit_wass_dist(punorm, d = 1)
    uwd1_spline <- unit_wass_dist(puspline, d = 1)
    uwd1_kern <- unit_wass_dist(pukern, d = 1)
    
    uwd2_cltn <- unit_wass_dist(pucltn, d = 2)
    uwd2_ordn <- unit_wass_dist(puordn, d = 2)
    uwd2_clt <- unit_wass_dist(puclt, d = 2)
    uwd2_ord <- unit_wass_dist(puord, d = 2)
    uwd2_ind <- unit_wass_dist(puind, d = 2)
    uwd2_meta <- unit_wass_dist(pumeta, d = 2)
    #uwd2_norm <- unit_wass_dist(punorm, d = 2)
    uwd2_spline <- unit_wass_dist(puspline, d = 2)
    uwd2_kern <- unit_wass_dist(pukern, d = 2)
    
    uwd1s <- c(uwd1_cltn, uwd1_ordn, uwd1_clt, uwd1_ord, uwd1_ind, uwd1_spline,
               uwd1_kern, uwd1_meta)
    uwd2s <- c(uwd2_cltn, uwd2_ordn, uwd2_clt, uwd2_ord, uwd2_ind, uwd2_spline,
               uwd2_kern, uwd2_meta)
    
     
    scores <- data.frame(rep = replicate, n = n, probs = p, quants = length(probs), 
               model = models, uwd1 = uwd1s, uwd2 = uwd2s)

    scores
  #write.csv(scores, "test_scores.csv", row.names = FALSE)  
    
    
  }


write.csv(distance, paste0("sim_scores/", dist, "5/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)



