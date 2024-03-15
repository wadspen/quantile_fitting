library(dplyr)
library(ggplot2)
library(cmdstanr)
library(rstanarm)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(doMC)

n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)


mod <- cmdstan_model(stan_file = 
                       '~/fitting_quantiles/fit_simple_normal_quantiles.stan')

ns <- c(50, 100, 200, 500, 1000, 5000)
mus <- c(-5, -2, 0, 2, 5)
sigmas <- c(.01, .1, .5, 1, 5, 10)

parameter_grid <- expand.grid(n = ns, mu = mus, sigma = sigmas)


fits <- foreach(i = 1:nrow(parameter_grid)) %dopar% {
    
    n <- parameter_grid$n[i]
    mu <- parameter_grid$mu[i]
    sigma <- parameter_grid$sigma[i]
    
    y <- rnorm(n, mu, sigma)
    probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    quantiles <- quantile(y, probs = probs, type = 3)
    
    dat <- list(
      N = length(quantiles),
      Q = quantiles,
      p = probs,
      n_components = 2,
      m = 0,
      c = 7,
      sv = 6,
      nv = 4000
    )
    
    samps <- mod$sample(data = dat,
                        iter_warmup = 500,
                        iter_sampling = 500,
                        chains = 1)
    
    diagnostics <- samps$diagnostic_summary()
    summary <- samps$summary()
    draws <- samps$draws(format = 'df')
    fit_list <- list(n = n, mu = mu, sigma = sigma,
                     diagnostics = diagnostics,
                     summary = summary,
                     draws = draws)
    
    fit_list
    
  }    


saveRDS(fits, "test_fits.rds")


















