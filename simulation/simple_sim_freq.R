#.libPaths("~/rlibs")
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
library(stringr)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
n.cores <- 62
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

args <- commandArgs()
dist <- args[6]
p <- as.numeric(args[7])
nind <- as.numeric(args[8])

print(dist); print(p); print(nind); print(model)


mod_loc <- "../stan_models/"



samples <- 50000
out_s <- 5000

samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
levels <- list(
  c(.25, .5, .75),
  c(.025, .25, .5, .95, .975),  
  seq(.1, .9, by = .1),
  c(.05, seq(.1, .9, by = .1), .95),
  c(.025, .05, seq(.1, .9, by = .1), .95, .975),
  c(.01, .025, .05, seq(.1, .9, by = .1), .95, .975, .99),
  seq(.05, .95, by = .05),
  c(.025, seq(.05, .95, by = .05), .975),
  c(.01, .025, seq(.05, .95, by = .05), .975, .99),
  seq(.01, .99, by = .02)	       
)



reps <- 200
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
                                  "VGAM", "distr", "dplyr", "LaplacesDemon",
                                  "stringr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
                      #foreach(n = samp_sizes, .combine = rbind) %:%
                      #foreach(p = 1:length(levels), .combine = rbind) %dopar% {
                      
                      set.seed(replicate)
                      source("./simulation_functions.R")
                      
                      if (dist == "norm") {
                        samp <- rnorm(n)
                        pdist <- function(x) {pnorm(x)}
                        ddist <- function(x) {dnorm(x)}
                        rdist <- function(n) {rnorm(n)}
                        qdist <- function(p) {qnorm(p)}
                      } else if (dist == "evd") {
                        samp <- revd(n)
                        pdist <- function(x) {pevd(x)}
                        ddist <- function(x) {devd(x)}
                        rdist <- function(n) {revd(n)}
                        qdist <- function(p) {qevd(p)}
                      } else if (dist == "lp") {
                        samp <- rlaplace(n)
                        pdist <- function(x) {plaplace(x)}
                        ddist <- function(x) {dlaplace(x)}
                        rdist <- function(n) {rlaplace(n)}
                        qdist <- function(p) {qlaplace(p)}
                      } else if (dist == "gmix") {
                        pars <- data.frame(mu = c(-1, 1.2),
                                           sigma = c(.9, .6),
                                           weight = c(.35, .65))
                        mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
                        samp <- r(mdist)(n)
                        pdist <- function(x) {p(mdist)(x)}
                        ddist <- function(x) {d(mdist)(x)}
                        rdist <- function(n) {r(mdist)(n)}
                        qdist <- function(p) {q(mdist)(p)}
                      } else if (dist == "tdp") {
                        Ndp <- 5
                        alphadp <- 1
                        
                        v <- rbeta(Ndp - 1, 1, alphadp)
                        pi <- Stick(v)
                        
                        pi <- c(pi, 1 - sum(pi))
                        pi[pi < 0] <- 0
                        while (sum(pi) != 1) {
                          pi <- pi/sum(pi)
                        }
                        
                        mus <- rnorm(Ndp + 1, 0, 2)
                        sigmas <- rep(1, Ndp + 1)
                        sigmas <- rgamma(Ndp + 1, 2)
                        
                        mdist <- make_normmix(mus, sigmas, pi)
                        samp <- r(mdist)(n)
                        pdist <- function(x) {p(mdist)(x)}
                        ddist <- function(x) {d(mdist)(x)}
                        rdist <- function(n) {r(mdist)(n)}
                        qdist <- function(p) {q(mdist)(p)}
                        
                      }
                      
                      probs <- levels[[p]]
                      quantiles <- quantile(samp, probs, type = 2)
                      true_quantiles <- qdist(probs)
                      
                      data <- data.frame(quantile = quantiles, prob = probs,
                                         true_quantile = true_quantiles)
                      
                      
                      start <- Sys.time() 
                      qspline <- make_q_fn(probs, quantiles)
                      rspline <- make_r_fn(probs, quantiles)
                      dspline <- make_d_fn(probs, quantiles)
                      puspline <- function(x) {pdist(qspline(x))}
                      end <- Sys.time()
                      tdiff <- end - start
                      tdiff_spline <- as.numeric(tdiff, "mins")
                      
                      
                      start <- Sys.time() 
                      rkern <- function(n) {rkden(n, quantiles, kernel = "gaussian")}
                      dkern <- function(x) {dkden(x, quantiles, kernel = "gaussian")}
                      
                      qkern <- function(p) {qkden(p, quantiles, kernel = "gaussian")}
                      pukern <- function(x) {pdist(qkern(x))}
                      end <- Sys.time()
                      tdiff <- end - start
                      tdiff_kern <- as.numeric(tdiff, "mins")
                      
                      
                      
                      
                      
                      
                      uwd1_spline <- unit_wass_dist(puspline, d = 1)
                      uwd1_kern <- unit_wass_dist(pukern, d = 1)
                      
                      ks_spline <- ks.test(pdist(rspline(out_s)), "punif")$statistic
                      ks_kern <- ks.test(rkern(out_s), "punif")$statistic
                      
                      kls <- rdist(samples)
                      splinex <- dspline(kls)
                      kernx <- dkern(kls)
                      py <- ddist(kls)
                      
                      spline_kl <- mean(log(py) - log(splinex))
                      kern_kl <- mean(log(py) - log(kernx))
                      
                      spline_tv <- dens_dist(dspline, ddist)
                      kern_tv <- dens_dist(dkern, ddist)
                     
          
                      
                      qdiff_spline <- sum(abs(qspline(probs) - qdist(probs))^2)
                      qdiff_kern <- sum(abs(qkern(probs) - qdist(probs))^2)
                      
                      data.frame(model = c("spline", "kern"), n = n, 
                                 numq = length(probs),
                                 dist = dist,
                                 samp_type = NA,
                                 rep = replicate,
                                 uwd1 = c(uwd1_spline, uwd1_kern), 
                                 kld = c(spline_kl, kern_kl),
                                 tv = c(spline_tv, kern_tv),
                                 comp95 = NA, 
                                 comp98 = NA,
                                 comp99 = NA,
                                 compg1 = NA,
                                 coverm = NA,
                                 qdiff = c(qdiff_spline, qdiff_kern),
                                 ftime = c(tdiff_spline, tdiff_kern)
                      )
                    }


saveRDS(distance, paste0("./simple_res/", 
                         paste(dist, "freq", p, nind, sep = "_"), ".rds"))



