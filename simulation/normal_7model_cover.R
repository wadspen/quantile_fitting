source("./simulation_functions.R")
library(cmdstanr)
library(distfromq)
library(evmix)
library(janitor)
library(distr)
library(tidyr)
library(parallel)
library(doParallel)
library(doMC)
# n.cores <- detectCores()
n.cores <- 9
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

args <- commandArgs()
sample_type = args[6]
p <- as.numeric(args[7])
nind <- as.numeric(args[8])
print(p); print(nind)

cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_normal_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_normal_ind_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_quantiles.stan')

cltnmod <- cmdstan_model(stan_file = 
                           '../stan_models/cdf_normal_n_quantiles.stan')

ordnmod <- cmdstan_model(stan_file = 
                           '../stan_models/order_normal_n_quantiles.stan')
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


mu <- 4
sigma <- 3.5
qtrue <- function(p) {qnorm(p, mu, sigma)}
ddist <- function(x) {dnorm(x, mu, sigma)}
rdist <- function(n) {rnorm(n, mu, sigma)}
out_s <- 10000
samples <- 50000
models <- c("cltn", "ordn", "clt", "ord", "ind", "spline", "kern")

reps <- 500
set.seed(92)
n <- samp_sizes[nind]
seeds <- sample(999999999, reps)
pdist <- function(x) {pnorm(x, mu, sigma)} 
qdist <- function(p) {qnorm(p, mu, sigma)}

distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "dplyr", "tidyr",
                                  "janitor")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %:% {
                      foreach(n = samp_sizes, .combine = rbind) %:%
                      foreach(p = 1:length(levels), .combine = rbind) %dopar% {
                      
                      true_params <- data.frame(variable = c("mu", "sigma", "n"), 
                                                truth = c(mu, sigma, n))
                      
                      source("./simulation_functions.R")
                      set.seed(replicate)
                      samp <- rnorm(n, mu, sigma)
                      probs <- levels[[p]]
                      quantiles <- quantile(samp, probs)
                      
                      data <- data.frame(quantile = quantiles, prob = probs,
                                         true_quantile = qdist(probs))
                      
                    
                      stan_data <- make_stan_data(data, size = n)
                      
                      # fit models
                      # if (mod %in% c("cltn", "ordn", "clt", "ord", "ind")) {
                      # if (mod == "cltn") {
                      # fit <- stan_fit_draws(cltnmod, stan_data)
                      # } else if (mod == "ordn") {
                      # fit <- stan_fit_draws(ordnmod, stan_data)
                      # } else if (mod == "clt") {
                      # fit <- stan_fit_draws(cltmod, stan_data)
                      # } else if (mod == "ord") {
                      # fit <- stan_fit_draws(ordmod, stan_data)
                      # } else if (mod == "ind") {
                      # fit <- stan_fit_draws(indmod, stan_data)
                      # }

                      # sum <- fit$summary
                      # sum_eval <- eval_sum(sum, true_params)
                      
                      ##############################################################
                      start <- Sys.time()
                      fit_cltn <- stan_fit_draws(cltnmod, stan_data, sampler = sample_type)
                      end <- Sys.time()
                      cltn_time <- difftime(end, start, units = "min")[[1]]
                      
                      draws <- fit_cltn[[2]]$draws(format = "df") %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      data_sum_cltn <- data %>% 
                        mutate(low95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               low90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.05),
                               upp90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.95),
                               mean = apply(draws, MARGIN = 2, 
                                            FUN = mean),
                               model = "cltn") 
           
                      ###############################################################
                      start <- Sys.time()
                      fit_ordn <- stan_fit_draws(ordnmod, stan_data, sampler = sample_type)
                      end <- Sys.time()
                      ordn_time <- difftime(end, start, units = "min")[[1]]
                      
                      draws <- fit_ordn[[2]]$draws(format = "df") %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      data_sum_ordn <- data %>% 
                        mutate(low95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               low90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.05),
                               upp90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.95),
                               mean = apply(draws, MARGIN = 2, 
                                            FUN = mean),
                               model = "ordn") 
                      
                      
                      ##############################################################
                      start <- Sys.time()
                      fit_clt <- stan_fit_draws(cltmod, stan_data, sampler = sample_type)
                      end <- Sys.time()
                      clt_time <- difftime(end, start, units = "min")[[1]]
                      
                      draws <- fit_clt[[2]]$draws(format = "df") %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      data_sum_clt <- data %>% 
                        mutate(low95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               low90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.05),
                               upp90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.95),
                               mean = apply(draws, MARGIN = 2, 
                                            FUN = mean),
                        model = "clt") 
                      
                      ##############################################################
                      start <- Sys.time()
                      fit_ord <- stan_fit_draws(ordmod, stan_data, sampler = sample_type)
                      end <- Sys.time()
                      ord_time <- difftime(end, start, units = "min")[[1]]
                      
                      draws <- fit_ord[[2]]$draws(format = "df") %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      data_sum_ord <- data %>% 
                        mutate(low95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               low90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.05),
                               upp90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.95),
                               mean = apply(draws, MARGIN = 2, 
                                            FUN = mean),
                               model = "ord") 
                      
                      
                      ##############################################################
                      start <- Sys.time()
                      fit_ind <- stan_fit_draws(indmod,stan_data, sampler = sample_type)
                      end <- Sys.time()
                      ind_time <- difftime(end, start, units = "min")[[1]]
                      
                      
                      
                      draws <- fit_ind[[2]]$draws(format = "df") %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      data_sum_ind <- data %>% 
                        mutate(low95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               low90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.05),
                               upp90 = apply(draws, MARGIN = 2, 
                                             FUN = quantile, probs = 0.95),
                               mean = apply(draws, MARGIN = 2, 
                                            FUN = mean),
                               model = "ind") 
                      
                      
                      
                      sum_data <- rbind(data_sum_cltn, data_sum_ordn, data_sum_clt,
                            data_sum_ord, data_sum_ind) %>%
                        rowwise() %>% 
                        mutate(cover = between(true_quantile, low95, 
                                               upp95)) %>% 
                        # group_by(model) %>% 
                        # summarise(mcover = mean(cover)) %>% 
                        mutate(rep = replicate, n = n, probs = p,
                               quants = length(probs)) %>% 
                        dplyr::select(rep, n, probs, quants, model,
                                      prob, quantile, true_quantile, low95,
                                      upp95, low90, upp90)
                      
                      
                      sum_data 
                      
                      
                      
                    }


write.csv(distance, paste0("sim_scores/norm_", sample_type, "/", 
                           "size", n, "_probs", length(levels[[p]]), 
                           "_coverage.csv"), row.names = FALSE)





