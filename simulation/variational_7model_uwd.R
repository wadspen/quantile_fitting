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
library(LaplacesDemon)
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
                          '../stan_models/cdf_quantile_normal_mix4.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_ind_quantile_normal_mix4.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_mix4_quantiles.stan')

cltnmod <- cmdstan_model(stan_file =
                           '../stan_models/cdf_quantile_normal_n_mix4.stan')

ordnmod <- cmdstan_model(stan_file =
                           '../stan_models/order_normal_n_mix4_quantiles.stan')

# metamod <- cmdstan_model(stan_file = 
#                            "../stan_models/metanorm_quantiles.stan")

#normmod <- cmdstan_model(stan_file =
#			   "../stan_models/normal_quantiles.stan")

mod_loc <- "../stan_models/"

burn <- 6000
samples <- 7000
out_s <- samples
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

models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern")


reps <- 1000
n <- samp_sizes[nind]


#p <- 2
#n <- 500
#dist <- "norm"
#replicate <- 2
distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                  "VGAM", "distr", "dplyr", "LaplacesDemon")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
  #foreach(n = samp_sizes, .combine = rbind) %:%
  #foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
   
    source("./simulation_functions.R")
    
    if (dist == "norm") {
      samp <- rnorm(n)
      pdist <- function(x) {pnorm(x)}
      ddist <- function(x) {dnorm(x)}
    } else if (dist == "evd") {
      samp <- revd(n)
      pdist <- function(x) {pevd(x)}
      ddist <- function(x) {devd(x)}
    } else if (dist == "lp") {
      samp <- rlaplace(n)
      pdist <- function(x) {plaplace(x)}
      ddist <- function(x) {dlaplace(x)}
    } else if (dist == "gmix") {
      pars <- data.frame(mu = c(-1, 1.2),
                         sigma = c(.9, .6),
                         weight = c(.35, .65))
      mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
      samp <- r(mdist)(n)
      pdist <- function(x) {p(mdist)(x)}
      ddist <- function(x) {d(mdist)(x)}
    }
    
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs, type = 2)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = n, comps = 4)
    
    
    #fit models
    fit_cltn <- stan_fit_draws(cltnmod, stan_data, 
                                sampler = "variational", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
	# saveRDS(list(draws = fit_cltn[[2]],samp = samp), paste0("sim_draws5/", dist, "_cltn_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
   # print("gets here") 
    fit_ordn <- stan_fit_draws(ordnmod, stan_data, 
                                sampler = "MCMC", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 

   	# saveRDS(list(draws = fit_ordn[[2]],samp = samp), paste0("sim_draws5/", dist, "_ordn_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
    fit_clt <- stan_fit_draws(cltmod, stan_data,
                               sampler = "MCMC", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 

   	# saveRDS(list(draws = fit_clt[[2]],samp = samp), paste0("sim_draws5/", dist, "_clt_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
    fit_ord <- stan_fit_draws(ordmod, stan_data,
                                sampler = "MCMC", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
    
   	# saveRDS(list(draws = fit_ord[[2]],samp = samp), paste0("sim_draws5/", dist, "_ord_rep", replicate, "_size", n, "_probs", length(probs), ".rds")) 
    fit_ind <- stan_fit_draws(indmod,stan_data,
                                sampler = "MCMC", burn = burn, samp = samples,
                                refresh = 100, out_s = 5000)
   
	# saveRDS(list(draws = fit_ind[[2]],samp = samp), paste0("sim_draws5/", dist, "_ind_rep", replicate, "_size", n, "_probs", length(probs), ".rds"))
    # fit_meta <- stan_fit_draws(metamod,stan_data,
    #                            sampler = "variational", burn = burn, samp = samples,
    #                             refresh = 100, out_s = 5000)
	# saveRDS(list(draws = fit_meta[[2]],samp = samp), paste0("sim_draws5/", dist, "_meta_rep", replicate, "_size", n, "_probs", length(probs), ".rds"))
#    fit_norm <- stan_fit_draws(normmod,stan_data,
#			       sampler = "variational",
#			       elbo = 300, grad = 10,
#			       out_s = 2000)
    
    params <- fit_cltn[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_mus <- m_params[1:4]
    cltn_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    cltn_wts <- wts/sum(wts)
    
    cltn_mix <- UnivarMixingDistribution(Norm(cltn_mus[1], cltn_sigmas[1]),
                                    Norm(cltn_mus[2], cltn_sigmas[2]),
                                    Norm(cltn_mus[3], cltn_sigmas[3]),
                                    Norm(cltn_mus[4], cltn_sigmas[4]),
                                    mixCoeff = cltn_wts)
    
    dcltn <- function(x) {d(cltn_mix)(x)}
    
    
    params <- fit_ordn[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    ordn_mus <- m_params[1:4]
    ordn_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    ordn_wts <- wts/sum(wts)
    
    ordn_mix <- UnivarMixingDistribution(Norm(ordn_mus[1], ordn_sigmas[1]),
                                    Norm(ordn_mus[2], ordn_sigmas[2]),
                                    Norm(ordn_mus[3], ordn_sigmas[3]),
                                    Norm(ordn_mus[4], ordn_sigmas[4]),
                                    mixCoeff = ordn_wts)
    
    dordn <- function(x) {d(ordn_mix)(x)}
    
    
    
    
    
    params <- fit_clt[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_mus <- m_params[1:4]
    clt_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    clt_wts <- wts/sum(wts)
    
    clt_mix <- UnivarMixingDistribution(Norm(clt_mus[1], clt_sigmas[1]),
                                    Norm(clt_mus[2], clt_sigmas[2]),
                                    Norm(clt_mus[3], clt_sigmas[3]),
                                    Norm(clt_mus[4], clt_sigmas[4]),
                                    mixCoeff = clt_wts)
    
    dclt <- function(x) {d(clt_mix)(x)}
    
    
    params <- fit_ord[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    ord_mus <- m_params[1:4]
    ord_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    ord_wts <- wts/sum(wts)
    
    ord_mix <- UnivarMixingDistribution(Norm(ord_mus[1], ord_sigmas[1]),
                                    Norm(ord_mus[2], ord_sigmas[2]),
                                    Norm(ord_mus[3], ord_sigmas[3]),
                                    Norm(ord_mus[4], ord_sigmas[4]),
                                    mixCoeff = ord_wts)
    
    dord <- function(x) {d(ord_mix)(x)}
    
    
    
    
    params <- fit_ind[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    ind_mus <- m_params[1:4]
    ind_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    ind_wts <- wts/sum(wts)
    
    ind_mix <- UnivarMixingDistribution(Norm(ind_mus[1], ind_sigmas[1]),
                                    Norm(ind_mus[2], ind_sigmas[2]),
                                    Norm(ind_mus[3], ind_sigmas[3]),
                                    Norm(ind_mus[4], ind_sigmas[4]),
                                    mixCoeff = ind_wts)
    
    dind <- function(x) {d(ind_mix)(x)}
    
    


    # fit_cltn <- fit_cltn[[1]]
    # fit_ordn <- fit_ordn[[1]]
    # fit_clt <- fit_clt[[1]]
    # fit_ord <- fit_ord[[1]]
    # fit_ind <- fit_ind[[1]]
    # fit_meta <- fit_meta[[1]]



    #unit draws
    udraws_cltn <- pdist(fit_cltn[[1]]$draws)
    udraws_ordn <- pdist(fit_ordn[[1]]$draws)
    udraws_clt <- pdist(fit_clt[[1]]$draws)
    udraws_ord <- pdist(fit_ord[[1]]$draws)
    udraws_ind <- pdist(fit_ind[[1]]$draws)
    # udraws_meta <- pdist(fit_meta$draws)
    #udraws_norm <- pdist(fit_norm$draws)
    
    
    #make unit ecdfs
    pucltn <- function(x) {ecdf(udraws_cltn)(x)}
    puordn <- function(x) {ecdf(udraws_ordn)(x)}
    puclt <- function(x) {ecdf(udraws_clt)(x)}
    puord <- function(x) {ecdf(udraws_ord)(x)}
    puind <- function(x) {ecdf(udraws_ind)(x)}
    # pumeta <- function(x) {ecdf(udraws_meta)(x)}
    punorm <- function(x) {ecdf(udraws_norm)(x)}
    qspline <- make_q_fn(probs, quantiles)
    rspline <- make_r_fn(probs, quantiles)
    dspline <- make_d_fn(probs, quantiles)
    rkern <- function(n) {rkden(n, quantiles, kernel = "epanechnikov")}
    dkern <- function(x) {dkden(x, quantiles, kernel = "epanechnikov")}
    puspline <- function(x) {pdist(qspline(x))}
    qkern <- function(p) {qkden(p, quantiles, kernel = "epanechnikov")}
    pukern <- function(x) {pdist(qkern(x))}
    
    
    uwd1_cltn <- unit_wass_dist(pucltn, d = 1)
    uwd1_ordn <- unit_wass_dist(puordn, d = 1)
    uwd1_clt <- unit_wass_dist(puclt, d = 1)
    uwd1_ord <- unit_wass_dist(puord, d = 1)
    uwd1_ind <- unit_wass_dist(puind, d = 1)
    # uwd1_meta <- unit_wass_dist(pumeta, d = 1)
    #uwd1_norm <- unit_wass_dist(punorm, d = 1)
    uwd1_spline <- unit_wass_dist(puspline, d = 1)
    uwd1_kern <- unit_wass_dist(pukern, d = 1)
    
    uwd2_cltn <- unit_wass_dist(pucltn, d = 2)
    uwd2_ordn <- unit_wass_dist(puordn, d = 2)
    uwd2_clt <- unit_wass_dist(puclt, d = 2)
    uwd2_ord <- unit_wass_dist(puord, d = 2)
    uwd2_ind <- unit_wass_dist(puind, d = 2)
    # uwd2_meta <- unit_wass_dist(pumeta, d = 2)
    #uwd2_norm <- unit_wass_dist(punorm, d = 2)
    uwd2_spline <- unit_wass_dist(puspline, d = 2)
    uwd2_kern <- unit_wass_dist(pukern, d = 2)
    
    
    ks_cltn <- ks.test(udraws_cltn, "punif")$statistic
    ks_ordn <- ks.test(udraws_ordn, "punif")$statistic
    ks_clt <- ks.test(udraws_clt, "punif")$statistic
    ks_ord <- ks.test(udraws_clt, "punif")$statistic
    ks_ind <- ks.test(udraws_ind, "punif")$statistic
    # ks_meta <- ks.test(udraws_meta, "punif")$statistic
    ks_spline <- ks.test(pdist(rspline(out_s)), "punif")$statistic
    ks_kern <- ks.test(rkern(out_s), "punif")$statistic
    
    
    unifs <- runif(5000)
    cltnx <- dcltn(unifs)
    ordnx <- dordn(unifs)
    cltx <- dclt(unifs)
    ordx <- dord(unifs)
    indx <- dind(unifs)
    splinex <- dspline(unifs)
    kernx <- dkern(unifs)
    py <- ddist(unifs)
    
    cltn_kl <- KLD(py, cltnx)$intrinsic.discrepancy
    ordn_kl <- KLD(py, ordnx)$intrinsic.discrepancy
    clt_kl <- KLD(py, cltx)$intrinsic.discrepancy
    ord_kl <- KLD(py, ordx)$intrinsic.discrepancy
    ind_kl <- KLD(py, indx)$intrinsic.discrepancy
    spline_kl <- KLD(py, splinex)$intrinsic.discrepancy
    kern_kl <- KLD(py, kernx)$intrinsic.discrepancy
    
    
    
    
    
    
    uwd1s <- c(uwd1_cltn, uwd1_ordn, uwd1_clt, uwd1_ord, uwd1_ind, uwd1_spline,
               uwd1_kern)
    uwd2s <- c(uwd2_cltn, uwd2_ordn, uwd2_clt, uwd2_ord, uwd2_ind, uwd2_spline,
               uwd2_kern)
    
    kss <- c(ks_cltn, ks_ordn, ks_clt, ks_ord, ks_ind, ks_spline, ks_kern)
    
    klds <- c(cltn_kl, ordn_kl, clt_kl, ord_kl, ind_kl, spline_kl, kern_kl)
    
     
    scores <- data.frame(rep = replicate, n = n, probs = p, quants = length(probs), 
               model = models, uwd1 = uwd1s, uwd2 = uwd2s, ks = kss, kld = klds)

    scores
  #write.csv(scores, "test_scores.csv", row.names = FALSE)  
    
    
  }


write.csv(distance, paste0("sim_scores/", dist, "/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)



