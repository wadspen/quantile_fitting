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
                          '../stan_models/cdf_quantile_normal_mixK_sb.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_normal_mixK_quantiles_sb.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_ind_quantile_normal_mixK_sb.stan')




burn <- 6000
sample <- 7000
samp_sizes <- c(50, 150, 500, 1000)
n <- 200
probs <- c(.01, .025, seq(.05, .95, by = .05), .975, .99)

#################################################
#################Dist fit images#################
#################################################
all_dens <- data.frame()
all_quants <- data.frame()
dist <- "lp"
for (n in samp_sizes) {

  if (dist == "lp") {
    samp <- rlaplace(n)
    true_quantiles <- qlaplace(probs)
  } else if (dist == "evd") {
    samp <- revd(n); dist <- "evd"
    true_quantiles <- qevd(probs)
  } else if (dist == "gmix") {
    pars <- data.frame(mu = c(-1, 1.2),
                       sigma = c(.9, .6),
                       weight = c(.35, .65))
    mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
    samp <- r(mdist)(n);
    true_quantiles <- q(mdist)(probs)
    ddist <- function(x) {d(mdist)(x)}
  }
  quantiles <- quantile(samp, probs = probs)
  
  
  dat <- data.frame(quantile = quantiles, prob = probs)
  stan_data <- make_stan_data(dat, size = n, comps = 12)
  
  
  cdfsamps <- cdfmod$sample(data = stan_data,
                                iter_warmup = burn,
                                iter_sampling = sample,
                                chains = 1,
                                # adapt_delta = .9999,
                                refresh = 1000)
  
  ordsamps <- ordmod$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            # adapt_delta = .9999,
                            refresh = 1000)
  
  
  indsamps <- indmod$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            # adapt_delta = .9999,
                            refresh = 1000)
  
  print("start processing data")
  print(paste("sample size is", n))
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

all_dens <- all_dens %>% 
  filter(row_number() > 20020)

all_quants <- all_quants %>% 
  filter(row_number() > 8284)
# saveRDS(all_dens, "../make_images/lp_dens_sb.rds")
# saveRDS(all_quants, "../make_images/lp_quants_sb.rds")
# saveRDS(all_dens, "../make_images/gmix_dens.rds")
# saveRDS(all_quants, "../make_images/gmix_quants.rds")

# all_dens <- readRDS("../make_images/lp_dens.rds")
# all_quants <- readRDS("../make_images/lp_quants.rds")
dens_fit <- all_dens %>% 
  # filter(row_number() > 20020) %>% 
  filter(model != "ord") %>%
  filter(N %in% c(50, 150, 500, 1000)) %>% 
  mutate(model = ifelse(model == "cdf", "QGP",
                        ifelse(model == "ind", "IND", 
                               ifelse(model == "kern", "KDE", 
                                      ifelse(model == "spline", "SPL",
                                             ifelse(model == "ord", "ORD", 
                                                    model)))))) %>% 
  mutate(model = factor(model, levels = c("KDE", "SPL", "IND", "QGP", "ORD"))) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = `0.025`, ymax = `0.975`),
              fill = "pink") + 
  geom_ribbon(aes(x = x, ymin = `0.25`, ymax = `0.75`),
              fill = "red") +
  geom_line(aes(x = x, y = yhat), colour = "darkgrey") +
  geom_line(aes(x = x, y = y)) + 
  ggh4x::facet_grid2(N~model,
                     scales = "free_y", independent = "y",
                     remove_labels = TRUE, ) +
  # facet_grid(N~model, scales = "free", independent = "all") +
  ylab("f(x)") +
  xlab("x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 11))


quant_fit <- all_quants %>% 
  # filter(row_number() > 8284) %>%
  filter(model != "ord") %>% 
  filter(N %in% c(50, 150, 500, 1000)) %>% 
  mutate(model = ifelse(model == "cdf", "QGP",
                        ifelse(model == "ind", "IND", 
                               ifelse(model == "kern", "KDE", 
                                      ifelse(model == "spline", "SPL",
                                             model))))) %>% 
  mutate(model = factor(model, levels = c("KDE", "SPL", "IND", "QGP"))) %>% 
  ggplot() +
  geom_ribbon(aes(x = prob, ymin = `0.025`, ymax = `0.975`),
              fill = "pink") + 
  geom_ribbon(aes(x = prob, ymin = `0.25`, ymax = `0.75`),
              fill = "red") +
  geom_line(aes(x = p, y = xhat), colour = "darkgrey", size = 1) +
  # geom_point(aes(x = prob, y = quantile)) +
  # geom_point(data = all_quants %>% 
  #              filter(!is.na(quantile)), aes(x = prob, y = quantile)) +
  # geom_point(data = data.frame(quantile = true_quantiles, prob = probs),
  #            aes(x = prob, y = quantile)) +
  geom_point(data = all_quants %>% 
               # filter(row_number() > 8284) %>%
               filter(!is.na(quantile), N <= 1000) %>% 
               dplyr::select(prob, quantile, N) %>% 
               unique(),
             aes(x = prob, y = quantile), size = .8) +
  ylab("Q(p)") +
  xlab("p") +
  # coord_cartesian(ylim = c(-5, 5)) + #laplace
  # coord_cartesian(ylim = c(-2.5, 4.2)) + #evd
  coord_cartesian(ylim = c(-3, 3.2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  ggh4x::facet_grid2(N~model,
                     # scales = "free", independent = "all",
                     remove_labels = TRUE, ) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 11))


cowplot::plot_grid(quant_fit, dens_fit, nrow = 1,
                   align = "hv")





















