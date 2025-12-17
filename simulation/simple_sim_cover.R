.libPaths("~/rlibs")
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
model <- args[9]

print(dist); print(p); print(nind); print(model)


mod_loc <- "../stan_models/"

mod_name <- ifelse(str_detect(model, "_shs"),
                   "cdf_quantile_normal_mix6_shs.stan",
                   "cdf_quantile_normal_mix6.stan")

mod <- cmdstan_model(stan_file = paste0(mod_loc, mod_name))





burn <- 2000
samples <- 2000
out_s <- 5000
sample_type <- "MCMC"

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
                                  "VGAM", "distr", "dplyr", "LaplacesDemon",
                                  "stringr")
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
                      }
                      
                      probs <- levels[[p]]
                      quantiles <- quantile(samp, probs, type = 2)
                      true_quantiles <- qdist(probs)
                      
                      data <- data.frame(quantile = quantiles, prob = probs,
                                         true_quantile = true_quantiles)
                      stan_data6 <- make_stan_data(data, size = n, comps = 6)
                      
                      
                      fit_mod <- stan_fit_draws(mod, stan_data6,
                                                 sampler = sample_type, burn = burn, samp = samples,
                                                 refresh = 100, out_s = 5000)
                      
                      
                      
                      fit <- fit_mod[[2]]
                      
                      draws <- fit$draws(format = "df") 
                      comps <- draws$samp_comp
                      comp_props <- table(comps)/sum(table(comps))
                      comp_props <- comp_props[order(comp_props, decreasing = TRUE)]
                      num_comps95 <- min(which(cumsum(comp_props) >= .95))
                      num_comps99 <- min(which(cumsum(comp_props) >= .99))
                      
                      
                      
                      drawsQ <- draws %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      
                      data_sum <- data %>% 
                        mutate(low95 = apply(drawsQ, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(drawsQ, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               mean = apply(drawsQ, MARGIN = 2, 
                                            FUN = mean)
                        ) 
                      
                      
                      coverm <- data_sum %>% 
                        mutate(cover = between(true_quantile, 
                                               low95, upp95)) %>% 
                        ungroup() %>% 
                        summarise(m = mean(cover)) %>% 
                        as.numeric()
                      
                      
                      
                      params <- draws %>% 
                        dplyr::select(contains("mus") | contains("sigmas") |
                                        contains("pi"))
                      m_params <- apply(params, MARGIN = 2, FUN = mean)
                      cltn_mus <- m_params[1:6]
                      cltn_sigmas <- m_params[7:12]
                      wts <- m_params[13:18] 
                      wts[wts < 0] <- 0
                      cltn_wts <- wts/sum(wts)
                      
                      
                      
                      
                      
                      
                      udraws <- pdist(fit_mod[[1]]$draws)
                      
                      
                      #make unit ecdfs
                      pu <- function(x) {ecdf(udraws)(x)}
                      
                      
                      
                      uwd1 <- unit_wass_dist(pu, d = 1)
                      
                      data.frame(model = model, n = n, 
                                 numq = length(probs),
                                 rep = replicate,
                                 uwd1 = uwd1, 
                                 comp95 = num_comps95, 
                                 comp99 = num_comps99,
                                 coverm = coverm
                                 )
                      # one <- uwd1_clt6
                    }


saveRDS(distance, paste0("./", model, ".rds"))
# write.csv(distance, paste0("sim_scores/", dist, "_", sample_type, "_numcomp",
#                            "/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)



