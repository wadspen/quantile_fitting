source("../simulation/simulation_functions.R")
source("../make_images/make_images_functions.R")
library(VGAM) #Laplace distribution
library(EnvStats) #Extreme value distribution
library(cmdstanr)
library(stringr)
library(ggplot2)
library(dplyr)
library(distr)
library(orderstats)
library(MASS) #for the multivariate normal model
library(tidyr)
library(evmix) #for the KDE estimation
library(distfromq)
library(ggpubr)
library(orderstats)



cdfmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_quantile_normal_mix4.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_mix4_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_ind_quantile_normal_mix4.stan')




burn <- 6000
sample <- 7000
samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
n <- 200
probs <- c(.01, .025, seq(.05, .95, by = .05), .975, .99)

#################################################
#################Dist fit images#################
#################################################
all_dens <- data.frame()
all_quants <- data.frame()
dist <- "lp"
for (n in samp_sizes) {

  samp <- rlaplace(n)
  true_quantiles <- qlaplace(probs)
  # samp <- revd(n)
  # true_quantiles <- qevd(probs)
  # pars <- data.frame(mu = c(-1, 1.2),
  #                    sigma = c(.9, .6),
  #                    weight = c(.35, .65))
  # mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
  # samp <- r(mdist)(n)
  # true_quantiles <- q(mdist)(probs)
  ddist <- function(x) {d(mdist)(x)}
  quantiles <- quantile(samp, probs = probs)
  
  
  dat <- data.frame(quantile = quantiles, prob = probs)
  stan_data <- make_stan_data(dat, size = n, comps = 4)
  
  
  cdfsamps <- cdfmod$sample(data = stan_data,
                                iter_warmup = burn,
                                iter_sampling = sample,
                                chains = 1,
                                # adapt_delta = .9999,
                                refresh = 100)
  
  ordsamps <- ordmod$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            # adapt_delta = .9999,
                            refresh = 100)
  
  
  indsamps <- indmod$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            # adapt_delta = .9999,
                            refresh = 100)
  
  print("start processing data")
  cdf_data <- get_quantile_samps(cdfsamps, quantiles = quantiles,
                                 n_modeled = TRUE,
                                 ind = FALSE, true_dist = dist, model = "cdf")
  
  print("cdf done")
  ord_data <- get_quantile_samps(ordsamps, quantiles = quantiles,
                                 n_modeled = TRUE, order = TRUE,
                                 ind = FALSE, 
                                 true_dist = dist, model = "ord")
  
  
  ind_data <- get_quantile_samps(indsamps, quantiles = quantiles, 
                                 true_dist = dist, model = "ind")
  
  print("data processed")
  dens <- rbind(cdf_data[[1]], ord_data[[1]], ind_data[[1]]) %>% 
    mutate(N = n)
  
  dkern <- function(x) {dkden(x, dat$quantile, kernel = "gaussian")}
  dspline <- make_d_fn(dat$prob, dat$quantile)
  
  # dx <- seq(-3.5, 3.5, length.out = 1001) #Laplace bounds
  # dx <- seq(-3, 8, length.out = 1001) #EVD bounds
  dx <- unique(dens$x) #GMIX
  dy <- dens$y[1:length(dx)]
  np_dens <- rbind(data.frame(x = dx, y = dy ,
                              yhat = dkern(dx), model = "kern", N = n),
        data.frame(x = dx, y = dy,
                   yhat = dspline(dx), model = "spline", N = n))
  
  dens <- dens %>% 
    full_join(np_dens, by = c("x", "y", "model", "N"))
  
  all_dens <- rbind(all_dens, dens)
  
  
  quants <- rbind(cdf_data[[2]], ord_data[[2]], ind_data[[2]]) %>% 
    mutate(N = n)
  
  
  qkern <- function(p) {qkden(p, dat$quantile, kernel = "gaussian")}
  qspline <- make_q_fn(dat$prob, dat$quantile)
  
  
  
  dp <- seq(.00001, .99999, length.out = 1001)
  np_quants <- rbind(data.frame(p = dp, xhat = qkern(dp),
                                model = "kern", N = n),
                     data.frame(p = dp, xhat = qspline(dp),
                                model = "spline", N = n))
  
  
  quants <- quants %>% 
    full_join(np_quants, by = c("model", "N"))
  
  all_quants <- rbind(all_quants, quants)

}

# saveRDS(all_dens, "./make_images/lp_dens.rds")
# saveRDS(all_quants, "./make_images/lp_quants.rds")

all_dens <- readRDS("../make_images/lp_dens.rds")
all_quants <- readRDS("../make_images/lp_quants.rds")
all_dens %>% 
  filter(model != "ord") %>% 
  filter(N %in% c(50, 150, 500, 1000)) %>% 
  ggplot() +
  geom_ribbon(aes(x = x, ymin = `0.025`, ymax = `0.975`),
              fill = "pink") + 
  geom_ribbon(aes(x = x, ymin = `0.25`, ymax = `0.75`),
              fill = "red") +
  geom_line(aes(x = x, y = yhat), colour = "darkgrey") +
  geom_line(aes(x = x, y = y)) + 
  ggh4x::facet_grid2(N~model,
                     scales = "free", independent = "all",
                     remove_labels = TRUE, ) +
  # facet_grid(N~model, scales = "free", independent = "all") +
  theme_bw() +
  theme(axis.text = element_blank())


all_quants %>% 
  filter(model != "ord") %>% 
  filter(N %in% c(50, 150, 500, 1000)) %>% 
  ggplot() +
  geom_ribbon(aes(x = prob, ymin = `0.025`, ymax = `0.975`),
              fill = "pink") + 
  geom_ribbon(aes(x = prob, ymin = `0.25`, ymax = `0.75`),
              fill = "red") +
  geom_line(aes(x = p, y = xhat), colour = "darkgrey") +
  geom_point(aes(x = prob, y = quantile)) +
  # geom_point(data = data.frame(quantile = true_quantiles, prob = probs),
  #            aes(x = prob, y = quantile)) +
  coord_cartesian(ylim = c(-5, 5)) +
  ggh4x::facet_grid2(N~model,
                     scales = "free", independent = "all",
                     remove_labels = TRUE, ) +
  theme_bw()






















