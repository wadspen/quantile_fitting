library(cmdstanr)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(stringr)
library(distr)
library(VGAM)

mod <- cmdstan_model(stan_file = 
                       './stan_models/fit_gaussian_mix5_quantiles.stan')

comp <- 5
n <- 400
probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

yln <- rlnorm(n)
ylp <- rlaplace(n)
yun <- runif(n)

quantilesln <- quantile(yln, probs = probs)
quantileslp <- quantile(ylp, probs = probs)
quantilesun <- quantile(yun, probs = probs)

datln <- list(
  N = length(quantilesln),
  Q = quantilesln,
  p = probs,
  n_components = 5,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

datlp <- list(
  N = length(quantileslp),
  Q = quantileslp,
  p = probs,
  n_components = 5,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

datun <- list(
  N = length(quantilesun),
  Q = quantilesun,
  p = probs,
  n_components = 5,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

sampsln <- mod$sample(data = datln,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)

sampslp <- mod$sample(data = datlp,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)




























