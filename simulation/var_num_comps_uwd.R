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
#n.cores <- detectCores()
n.cores <- 64
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


cltmod5 <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_quantile_normal_mix5.stan')


cltnmod5 <- cmdstan_model(stan_file =
                           '../stan_models/cdf_quantile_normal_n_mix5.stan')


cltmod4 <- cmdstan_model(stan_file = 
                           '../stan_models/cdf_quantile_normal_mix4.stan')


cltnmod4 <- cmdstan_model(stan_file =
                            '../stan_models/cdf_quantile_normal_n_mix4.stan')


cltmod3 <- cmdstan_model(stan_file = 
                           '../stan_models/cdf_quantile_normal_mix3.stan')


cltnmod3 <- cmdstan_model(stan_file =
                            '../stan_models/cdf_quantile_normal_n_mix3.stan')

cltmod2 <- cmdstan_model(stan_file = 
                           '../stan_models/cdf_quantile_normal_mix2.stan')


cltnmod2 <- cmdstan_model(stan_file =
                            '../stan_models/cdf_quantile_normal_n_mix2.stan')


cltmod1 <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_normal_quantiles.stan')


cltnmod1 <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_normal_n_quantiles.stan')





#normmod <- cmdstan_model(stan_file =
#			   "../stan_models/normal_quantiles.stan")

mod_loc <- "../stan_models/"

burn <- 20000
samples <- 60000
out_s <- 5000
sample_type <- "MCMC"
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

models <- c("cltn5", "clt5", "cltn4", "clt4", "cltn3", 
            "clt3", "cltn2", "clt2", "cltn1", "clt1")

comps <- rep(5:1, each = 2)


reps <- 500
n <- samp_sizes[nind]

#p <- 7
#n <- 1000
#replicate <- 4
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
      rdist <- function(n) {rnorm(n)}
    } else if (dist == "evd") {
      samp <- revd(n)
      pdist <- function(x) {pevd(x)}
      ddist <- function(x) {devd(x)}
      rdist <- function(n) {revd(n)}
    } else if (dist == "lp") {
      samp <- rlaplace(n)
      pdist <- function(x) {plaplace(x)}
      ddist <- function(x) {dlaplace(x)}
      rdist <- function(n) {rlaplace(n)}
    } else if (dist == "gmix") {
      pars <- data.frame(mu = c(-1, 1.2),
                         sigma = c(.9, .6),
                         weight = c(.35, .65))
      mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
      samp <- r(mdist)(n)
      pdist <- function(x) {p(mdist)(x)}
      ddist <- function(x) {d(mdist)(x)}
      rdist <- function(n) {r(mdist)(n)}
    }
    
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs, type = 2)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data5 <- make_stan_data(data, size = n, comps = 5)
    
    
    #fit models
    fit_cltn5 <- stan_fit_draws(cltnmod5, stan_data5, 
                                sampler = sample_type, burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 

    fit_clt5 <- stan_fit_draws(cltmod5, stan_data5,
                               sampler = sample_type, burn = burn, samp = samples,
                                refresh = 100, out_s = 5000)
    
    
    
    
    params <- fit_cltn5[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_mus <- m_params[1:5]
    cltn_sigmas <- m_params[6:10]
    wts <- m_params[11:15] 
    wts[wts < 0] <- 0
    cltn_wts <- wts/sum(wts)
    
    cltn_mix5 <- UnivarMixingDistribution(Norm(cltn_mus[1], cltn_sigmas[1]),
                                         Norm(cltn_mus[2], cltn_sigmas[2]),
                                         Norm(cltn_mus[3], cltn_sigmas[3]),
                                         Norm(cltn_mus[4], cltn_sigmas[4]),
                                         Norm(cltn_mus[5], cltn_sigmas[5]),
                                         mixCoeff = cltn_wts)
    
    dcltn5 <- function(x) {d(cltn_mix5)(x)}
    qcltn5 <- function(p) {q(cltn_mix5)(p)} 
    
    
    
    
    
    
    params <- fit_clt5[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_mus <- m_params[1:5]
    clt_sigmas <- m_params[6:10]
    wts <- m_params[11:15] 
    wts[wts < 0] <- 0
    clt_wts <- wts/sum(wts)
    
    clt_mix5 <- UnivarMixingDistribution(Norm(clt_mus[1], clt_sigmas[1]),
                                        Norm(clt_mus[2], clt_sigmas[2]),
                                        Norm(clt_mus[3], clt_sigmas[3]),
                                        Norm(clt_mus[4], clt_sigmas[4]),
                                        Norm(clt_mus[5], clt_sigmas[5]),
                                        mixCoeff = clt_wts)
    
    dclt5 <- function(x) {d(clt_mix5)(x)}
    qclt5 <- function(p) {q(clt_mix5)(p)} 
    
    
    ####4 components
    stan_data4 <- make_stan_data(data, size = n, comps = 4)
    
    
    #fit models
    fit_cltn4 <- stan_fit_draws(cltnmod4, stan_data4, 
                                sampler = sample_type, burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
    
    fit_clt4 <- stan_fit_draws(cltmod4, stan_data4,
                               sampler = sample_type, burn = burn, samp = samples,
                               refresh = 100, out_s = 5000)
    
    
    
    params <- fit_cltn4[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_mus <- m_params[1:4]
    cltn_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    cltn_wts <- wts/sum(wts)
    
    cltn_mix4 <- UnivarMixingDistribution(Norm(cltn_mus[1], cltn_sigmas[1]),
                                         Norm(cltn_mus[2], cltn_sigmas[2]),
                                         Norm(cltn_mus[3], cltn_sigmas[3]),
                                         Norm(cltn_mus[4], cltn_sigmas[4]),
                                         # Norm(cltn_mus[5], cltn_sigmas[5]),
                                         mixCoeff = cltn_wts)
    
    dcltn4 <- function(x) {d(cltn_mix4)(x)}
    qcltn4 <- function(p) {q(cltn_mix4)(p)} 
    
    
    
    
    
    
    params <- fit_clt4[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_mus <- m_params[1:4]
    clt_sigmas <- m_params[5:8]
    wts <- m_params[9:12] 
    wts[wts < 0] <- 0
    clt_wts <- wts/sum(wts)
    
    clt_mix4 <- UnivarMixingDistribution(Norm(clt_mus[1], clt_sigmas[1]),
                                        Norm(clt_mus[2], clt_sigmas[2]),
                                        Norm(clt_mus[3], clt_sigmas[3]),
                                        Norm(clt_mus[4], clt_sigmas[4]),
                                        # Norm(clt_mus[5], clt_sigmas[5]),
                                        mixCoeff = clt_wts)
    
    dclt4 <- function(x) {d(clt_mix4)(x)}
    qclt4 <- function(p) {q(clt_mix4)(p)} 
    
    ####3 components
    stan_data3 <- make_stan_data(data, size = n, comps = 3)
    
    
    #fit models
    fit_cltn3 <- stan_fit_draws(cltnmod3, stan_data3, 
                                sampler = sample_type, burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
    
    fit_clt3 <- stan_fit_draws(cltmod3, stan_data3,
                               sampler = sample_type, burn = burn, samp = samples,
                               refresh = 100, out_s = 5000)
    
    
    
    params <- fit_cltn3[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_mus <- m_params[1:3]
    cltn_sigmas <- m_params[4:6]
    wts <- m_params[7:9] 
    wts[wts < 0] <- 0
    cltn_wts <- wts/sum(wts)
    
    cltn_mix3 <- UnivarMixingDistribution(Norm(cltn_mus[1], cltn_sigmas[1]),
                                         Norm(cltn_mus[2], cltn_sigmas[2]),
                                         Norm(cltn_mus[3], cltn_sigmas[3]),
                                         # Norm(cltn_mus[4], cltn_sigmas[4]),
                                         # Norm(cltn_mus[5], cltn_sigmas[5]),
                                         mixCoeff = cltn_wts)
    
    dcltn3 <- function(x) {d(cltn_mix3)(x)}
    qcltn3 <- function(p) {q(cltn_mix3)(p)} 
    
    
    
    
    
    
    params <- fit_clt3[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_mus <- m_params[1:3]
    clt_sigmas <- m_params[4:6]
    wts <- m_params[7:9] 
    wts[wts < 0] <- 0
    clt_wts <- wts/sum(wts)
    
    clt_mix3 <- UnivarMixingDistribution(Norm(clt_mus[1], clt_sigmas[1]),
                                        Norm(clt_mus[2], clt_sigmas[2]),
                                        Norm(clt_mus[3], clt_sigmas[3]),
                                        # Norm(clt_mus[4], clt_sigmas[4]),
                                        # Norm(clt_mus[5], clt_sigmas[5]),
                                        mixCoeff = clt_wts)
    
    dclt3 <- function(x) {d(clt_mix3)(x)}
    qclt3 <- function(p) {q(clt_mix3)(p)} 
    
    
    
    ####2 components
    stan_data2 <- make_stan_data(data, size = n, comps = 2)
    
    
    #fit models
    fit_cltn2 <- stan_fit_draws(cltnmod2, stan_data2, 
                                sampler = sample_type, burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
    
    fit_clt2 <- stan_fit_draws(cltmod2, stan_data2,
                               sampler = sample_type, burn = burn, samp = samples,
                               refresh = 100, out_s = 5000)
    
    params <- fit_cltn2[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_mus <- m_params[1:2]
    cltn_sigmas <- m_params[3:4]
    wts <- m_params[5:6] 
    wts[wts < 0] <- 0
    cltn_wts <- wts/sum(wts)
    
    cltn_mix2 <- UnivarMixingDistribution(Norm(cltn_mus[1], cltn_sigmas[1]),
                                         Norm(cltn_mus[2], cltn_sigmas[2]),
                                         # Norm(cltn_mus[3], cltn_sigmas[3]),
                                         # Norm(cltn_mus[4], cltn_sigmas[4]),
                                         # Norm(cltn_mus[5], cltn_sigmas[5]),
                                         mixCoeff = cltn_wts)
    
    dcltn2 <- function(x) {d(cltn_mix2)(x)}
    qcltn2 <- function(p) {q(cltn_mix2)(p)} 
    
    
    
    
    
    
    params <- fit_clt2[[2]]$draws(variables = c("mus", "sigmas", "pi"), 
                                 format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_mus <- m_params[1:2]
    clt_sigmas <- m_params[3:4]
    wts <- m_params[5:6] 
    wts[wts < 0] <- 0
    clt_wts <- wts/sum(wts)
    
    clt_mix2 <- UnivarMixingDistribution(Norm(clt_mus[1], clt_sigmas[1]),
                                        Norm(clt_mus[2], clt_sigmas[2]),
                                        # Norm(clt_mus[3], clt_sigmas[3]),
                                        # Norm(clt_mus[4], clt_sigmas[4]),
                                        # Norm(clt_mus[5], clt_sigmas[5]),
                                        mixCoeff = clt_wts)
    
    dclt2 <- function(x) {d(clt_mix2)(x)}
    qclt2 <- function(p) {q(clt_mix2)(p)} 
    
    
    ####1 component
    stan_data1 <- make_stan_data(data, size = n, comps = 1)
    
    
    #fit models
    fit_cltn1 <- stan_fit_draws(cltnmod1, stan_data1, 
                                sampler = sample_type, burn = burn, samp = samples,
                                refresh = 100, out_s = 5000) 
    
    fit_clt1 <- stan_fit_draws(cltmod1, stan_data1,
                               sampler = sample_type, burn = burn, samp = samples,
                               refresh = 100, out_s = 5000)


    
    
    
    params <- fit_cltn1[[2]]$draws(variables = c("mu", "sigma"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    cltn_mus <- m_params[1]
    cltn_sigmas <- m_params[2]
    wts <- 1
    wts[wts < 0] <- 0
    cltn_wts <- wts/sum(wts)
    
    cltn_mix1 <- UnivarMixingDistribution(Norm(cltn_mus[1], cltn_sigmas[1]),
#                                     Norm(cltn_mus[2], cltn_sigmas[2]),
#                                     Norm(cltn_mus[3], cltn_sigmas[3]),
#                                     Norm(cltn_mus[4], cltn_sigmas[4]),
# 				    Norm(cltn_mus[5], cltn_sigmas[5]),
                                    mixCoeff = cltn_wts)
    
    dcltn1 <- function(x) {d(cltn_mix1)(x)}
    qcltn1 <- function(p) {q(cltn_mix1)(p)} 
    

    
    
    
    
    params <- fit_clt1[[2]]$draws(variables = c("mu", "sigma"), 
                                  format = "df")               
    m_params <- apply(params, MARGIN = 2, FUN = mean)
    clt_mus <- m_params[1]
    clt_sigmas <- m_params[2]
    wts <- 1
    wts[wts < 0] <- 0
    clt_wts <- wts/sum(wts)
    
    clt_mix1 <- UnivarMixingDistribution(Norm(clt_mus[1], clt_sigmas[1]),
#                                     Norm(clt_mus[2], clt_sigmas[2]),
#                                     Norm(clt_mus[3], clt_sigmas[3]),
#                                     Norm(clt_mus[4], clt_sigmas[4]),
# 				    Norm(clt_mus[5], clt_sigmas[5]),
                                    mixCoeff = clt_wts)
    
    dclt1 <- function(x) {d(clt_mix1)(x)}
    qclt1 <- function(p) {q(clt_mix1)(p)} 
    




    #unit draws
    udraws_cltn5 <- pdist(fit_cltn5[[1]]$draws)
    udraws_clt5 <- pdist(fit_clt5[[1]]$draws)
    
    udraws_cltn4 <- pdist(fit_cltn4[[1]]$draws)
    udraws_clt4 <- pdist(fit_clt4[[1]]$draws)
    
    udraws_cltn3 <- pdist(fit_cltn3[[1]]$draws)
    udraws_clt3 <- pdist(fit_clt3[[1]]$draws)
    
    udraws_cltn2 <- pdist(fit_cltn2[[1]]$draws)
    udraws_clt2 <- pdist(fit_clt2[[1]]$draws)
    
    udraws_cltn1 <- pdist(fit_cltn1[[1]]$draws)
    udraws_clt1 <- pdist(fit_clt1[[1]]$draws)
    
    
    #make unit ecdfs
    pucltn5 <- function(x) {ecdf(udraws_cltn5)(x)}
    puclt5 <- function(x) {ecdf(udraws_clt5)(x)}
    
    pucltn4 <- function(x) {ecdf(udraws_cltn4)(x)}
    puclt4 <- function(x) {ecdf(udraws_clt4)(x)}
    
    pucltn3 <- function(x) {ecdf(udraws_cltn3)(x)}
    puclt3 <- function(x) {ecdf(udraws_clt3)(x)}
    
    pucltn2 <- function(x) {ecdf(udraws_cltn2)(x)}
    puclt2 <- function(x) {ecdf(udraws_clt2)(x)}
    
    pucltn1 <- function(x) {ecdf(udraws_cltn1)(x)}
    puclt1 <- function(x) {ecdf(udraws_clt1)(x)}
    
    
    
    
    
    
    pucltno5 <- function(x) {pdist(qcltn5(x))}
    puclto5 <- function(x) {pdist(qclt5(x))}
    
    pucltno4 <- function(x) {pdist(qcltn4(x))}
    puclto4 <- function(x) {pdist(qclt4(x))}
    
    pucltno3 <- function(x) {pdist(qcltn3(x))}
    puclto3 <- function(x) {pdist(qclt3(x))}
    
    pucltno2 <- function(x) {pdist(qcltn2(x))}
    puclto2 <- function(x) {pdist(qclt2(x))}
    
    pucltno1 <- function(x) {pdist(qcltn1(x))}
    puclto1 <- function(x) {pdist(qclt1(x))}
    
    
    
    
    uwd1_cltn5 <- unit_wass_dist(pucltn5, d = 1)
    uwd1_clt5 <- unit_wass_dist(puclt5, d = 1)
    
    uwd1_cltn4 <- unit_wass_dist(pucltn4, d = 1)
    uwd1_clt4 <- unit_wass_dist(puclt4, d = 1)
    
    uwd1_cltn3 <- unit_wass_dist(pucltn3, d = 1)
    uwd1_clt3 <- unit_wass_dist(puclt3, d = 1)
    
    uwd1_cltn2 <- unit_wass_dist(pucltn2, d = 1)
    uwd1_clt2 <- unit_wass_dist(puclt2, d = 1)
    
    uwd1_cltn1 <- unit_wass_dist(pucltn1, d = 1)
    uwd1_clt1 <- unit_wass_dist(puclt1, d = 1)
    
    
    
    
    
    
    
    uwd1_cltno5 <- unit_wass_dist(pucltno5, d = 1)
    uwd1_clto5 <- unit_wass_dist(puclto5, d = 1)
    
    uwd1_cltno4 <- unit_wass_dist(pucltno4, d = 1)
    uwd1_clto4 <- unit_wass_dist(puclto4, d = 1)
    
    uwd1_cltno3 <- unit_wass_dist(pucltno3, d = 1)
    uwd1_clto3 <- unit_wass_dist(puclto3, d = 1)
    
    uwd1_cltno2 <- unit_wass_dist(pucltno2, d = 1)
    uwd1_clto2 <- unit_wass_dist(puclto2, d = 1)
    
    uwd1_cltno1 <- unit_wass_dist(pucltno1, d = 1)
    uwd1_clto1 <- unit_wass_dist(puclto1, d = 1)
   

    
    
    
    uwd2_cltn5 <- unit_wass_dist(pucltn5, d = 2)
    uwd2_clt5 <- unit_wass_dist(puclt5, d = 2)
    
    uwd2_cltn4 <- unit_wass_dist(pucltn4, d = 2)
    uwd2_clt4 <- unit_wass_dist(puclt4, d = 2)
    
    uwd2_cltn3 <- unit_wass_dist(pucltn3, d = 2)
    uwd2_clt3 <- unit_wass_dist(puclt3, d = 2)
    
    uwd2_cltn2 <- unit_wass_dist(pucltn2, d = 2)
    uwd2_clt2 <- unit_wass_dist(puclt2, d = 2)
    
    uwd2_cltn1 <- unit_wass_dist(pucltn1, d = 2)
    uwd2_clt1 <- unit_wass_dist(puclt1, d = 2)
    
    
    
    
    ks_cltn5 <- ks.test(udraws_cltn5, "punif")$statistic
    ks_clt5 <- ks.test(udraws_clt5, "punif")$statistic
    
    ks_cltn4 <- ks.test(udraws_cltn4, "punif")$statistic
    ks_clt4 <- ks.test(udraws_clt4, "punif")$statistic
    
    ks_cltn3 <- ks.test(udraws_cltn3, "punif")$statistic
    ks_clt3 <- ks.test(udraws_clt3, "punif")$statistic
    
    ks_cltn2 <- ks.test(udraws_cltn2, "punif")$statistic
    ks_clt2 <- ks.test(udraws_clt2, "punif")$statistic
    
    ks_cltn1 <- ks.test(udraws_cltn1, "punif")$statistic
    ks_clt1 <- ks.test(udraws_clt1, "punif")$statistic
    
    
    
    kls <- rdist(samples)
    cltnx5 <- dcltn5(kls)
    cltx5 <- dclt5(kls)
    
    cltnx4 <- dcltn4(kls)
    cltx4 <- dclt4(kls)
    
    cltnx3 <- dcltn3(kls)
    cltx3 <- dclt3(kls)
    
    cltnx2 <- dcltn2(kls)
    cltx2 <- dclt2(kls)
    
    cltnx1 <- dcltn1(kls)
    cltx1 <- dclt1(kls)
    
    py <- ddist(kls)
    
    cltn_kl5 <- mean(log(py) - log(cltnx5))
    clt_kl5 <- mean(log(py) - log(cltx5))
    
    cltn_kl4 <- mean(log(py) - log(cltnx4))
    clt_kl4 <- mean(log(py) - log(cltx4))
    
    cltn_kl3 <- mean(log(py) - log(cltnx3))
    clt_kl3 <- mean(log(py) - log(cltx3))
    
    cltn_kl2 <- mean(log(py) - log(cltnx2))
    clt_kl2 <- mean(log(py) - log(cltx2))
    
    cltn_kl1 <- mean(log(py) - log(cltnx1))
    clt_kl1 <- mean(log(py) - log(cltx1))
    
    
    
    
    cltn_tv5 <- dens_dist(dcltn5, ddist)
    clt_tv5 <- dens_dist(dclt5, ddist)
    
    cltn_tv4 <- dens_dist(dcltn4, ddist)
    clt_tv4 <- dens_dist(dclt4, ddist)
    
    cltn_tv3 <- dens_dist(dcltn3, ddist)
    clt_tv3 <- dens_dist(dclt3, ddist)
    
    cltn_tv2 <- dens_dist(dcltn2, ddist)
    clt_tv2 <- dens_dist(dclt2, ddist)
    
    cltn_tv1 <- dens_dist(dcltn1, ddist)
    clt_tv1 <- dens_dist(dclt1, ddist)
    
    
    
    
    uwd1s <- c(uwd1_cltn5, uwd1_clt5, uwd1_cltn4, uwd1_clt4, uwd1_cltn3, 
               uwd1_clt3,
               uwd1_cltn2, uwd1_clt2, uwd1_cltn1, uwd1_clt1)
    
    uwd2s <- c(uwd2_cltn5, uwd2_clt5, uwd2_cltn4, uwd2_clt4, 
               uwd2_cltn3, uwd2_clt3, uwd2_cltn2, uwd2_clt2,
               uwd2_cltn1, uwd2_clt1)

    uwd1os <- c(uwd1_cltno5, uwd1_clto5, uwd1_cltno4, uwd1_clto4,
                uwd1_cltno3, uwd1_clto3, uwd1_cltno2, uwd1_clto2,
                uwd1_cltno1, uwd1_clto1) 
    
    kss <- c(ks_cltn5, ks_clt5, ks_cltn4, ks_clt4, ks_cltn3, ks_clt3,
             ks_cltn2, ks_clt2, ks_cltn1, ks_clt1)
    
    klds <- c(cltn_kl5, clt_kl5, cltn_kl4, clt_kl4, cltn_kl3, clt_kl3,
              cltn_kl2, clt_kl2, cltn_kl1, clt_kl1)
    
    tvs <- c(cltn_tv5, clt_tv5, cltn_tv4, clt_tv4, cltn_tv3, clt_tv3,
             cltn_tv2, clt_tv2, cltn_tv1, clt_tv1)
    
     
    scores <- data.frame(rep = replicate, n = n, probs = p, quants = length(probs), 
               model = models, comp = comps, uwd1 = uwd1s, uwd1o = uwd1os, 
               uwd2 = uwd2s, ks = kss, kld = klds,
               tv = tvs)

    scores
  #write.csv(scores, "test_scores.csv", row.names = FALSE)  
    
    

}

#print(paste0("sim_scores/", dist, "_", sample_type, "/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)
#write.csv(distance, "test_distance.csv")
write.csv(distance, paste0("sim_scores/", dist, "_", sample_type, "_numcomp",
                           "/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)



