source("./simulation/simulation_functions.R")
source("make_images_functions.R")
library(VGAM) #Laplace distribution
library(EnvStats) #Extreme value distribution
library(cmdstanr)
library(stringr)
library(ggplot2)
library(dplyr)
library(distr)
library(MASS) #for the multivariate normal model
library(tidyr)
library(evmix) #for the KDE estimation
library(distfromq)
library(ggpubr)



cdfmod <- cmdstan_model(stan_file = 
                          './stan_models/cdf_quantile_normal_mix4.stan')

ordmod <- cmdstan_model(stan_file = 
                          './stan_models/order_normal_mix4_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          './stan_models/cdf_ind_quantile_normal_mix4.stan')




burn <- 600
sample <- 700

n <- 200
probs <- c(.01, .025, seq(.05, .95, by = .05), .975, .99)

#################################################
#################Dist fit images#################
#################################################
samp <- rlaplace(n)
true_quantiles <- qlaplace(probs)
# samp <- revd(n)
# pars <- data.frame(mu = c(-1, 1.2),
#                    sigma = c(.9, .6),
#                    weight = c(.35, .65))
# mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
# samp <- r(mdist)(n)
# true_quantiles <- q(mdist)(probs)
# ddist <- function(x) {d(mdist)(x)}
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


cdf_data <- get_quantile_samps(cdfsamps, quantiles = true_quantiles,
                               n_modeled = TRUE,
                               ind = FALSE, true_dist = dist)
ord_data <- get_quantile_samps(ordsamps, quantiles = true_quantiles,
                               n_modeled = TRUE, order = TRUE,
                               ind = FALSE, 
                               true_dist = dist)
ind_data <- get_quantile_samps(indsamps, quantiles = true_quantiles, 
                               true_dist = dist)







cdf_quant <- make_quant_plot(cdf_data)
ord_quant <- make_quant_plot(ord_data)
ind_quant <- make_quant_plot(ind_data)


qkern <- function(p) {qkden(p, dat$quantile, kernel = "epanechnikov")}
qspline <- make_q_fn(dat$prob, dat$quantile)

qp <- seq(.0001, .9999, length.out = 1001)
est_quant_np <- data.frame(p = qp, qkerny = qkern(qp), 
           qspliney <- qspline(qp))

quant <- dat %>% 
  ggplot() +
  geom_point(aes(x = prob, y = quantile), size = 1.3) +
  xlab("") +
  ylab("") +
  # coord_cartesian(ylim=c(-2.6, 5.5)) + #EVD
  # coord_cartesian(ylim=c(-5, 5)) + #Laplace
  coord_cartesian(ylim=c(-2.5, 3)) + #gmix
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=23),
        plot.title=element_text(hjust = .1, margin = margin(t = 5, b=-20)))


kern_quant <- dat %>% 
  ggplot() +
  geom_line(data = est_quant_np, aes(x = p, y = qkerny), size = .75,
            colour = "darkgrey") +
  geom_point(aes(x = prob, y = quantile), size = 1.3) +
  xlab("") +
  ylab("") +
  # coord_cartesian(ylim=c(-2.6, 5.5)) + #EVD
  # coord_cartesian(ylim=c(-5, 5)) + #Laplace
  coord_cartesian(ylim=c(-2.5, 3)) + #gmix
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=23),
        plot.title=element_text(hjust = .1, margin = margin(t = 5, b=-20)))


spline_quant <- dat %>%
  ggplot() +
  geom_line(data = est_quant_np, aes(x = p, y = qspliney), size = .75,
            colour = "darkgrey") +
  geom_point(aes(x = prob, y = quantile), size = 1.3) +
  xlab("") +
  ylab("") +
  # coord_cartesian(ylim=c(-2.6, 5.5)) + #EVD
  # coord_cartesian(ylim=c(-5, 5)) + #Laplace
  coord_cartesian(ylim=c(-2.5, 3)) + #gmix
  # ggtitle("(c)") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=23),
        plot.title=element_text(hjust = .1, margin = margin(t = 0, b=-20)))

margin = theme(plot.margin = unit(c(.1,.2,.1,.1), "cm"))
as_ggplot(gridExtra::grid.arrange(grobs = 
                                    lapply(list(quant + ggtitle("(a)"), 
                                  kern_quant + ggtitle("(b)"),
                                  spline_quant + ggtitle("(c)"), 
                                  ind_quant + ggtitle("(d)"),
                                  ord_quant + ggtitle("(e)"), 
                                  cdf_quant + ggtitle("(f)")), "+", margin), 
                                  left = text_grob("Q(p)", size = 23, 
                                                   vjust = 1, rot = 90),
                                  bottom = text_grob("p", size = 23,
                                                     vjust = -.5)))






#####################################################################
##########################Density plots##############################
#####################################################################


cdf_dens <- make_dens_plot(cdf_data)
ord_dens <- make_dens_plot(ord_data)
ind_dens <- make_dens_plot(ind_data)


dkern <- function(x) {dkden(x, dat$quantile, kernel = "epanechnikov")}
dspline <- make_d_fn(dat$prob, dat$quantile)

# dx <- seq(-3.5, 3.5, length.out = 1001) #Laplace bounds
# dx <- seq(-3, 8, length.out = 1001) #EVD bounds
dx <- seq(-3, 3, length.out = 1001) #GMIX
est_dens_np <- data.frame(x = dx, dkerny = dkern(dx), 
                           dspliney = dspline(dx),
                           # y = dlaplace(dx)
                           # y = devd(dx)
                           y = ddist(dx)
                          )



dens <- est_dens_np %>% 
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=23),
        plot.title=element_text(hjust = .1, margin = margin(t = 5, b=-20)))


kern_dens <- est_dens_np %>% 
  ggplot() +
  geom_line(aes(x = x, y = dkerny), size = 1.3, colour = "darkgrey") +
  geom_line(aes(x = x, y = y), size = .8) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=23),
        plot.title=element_text(hjust = .1, margin = margin(t = 5, b=-20)))


spline_dens <- est_dens_np %>% 
  ggplot() +
  geom_line(aes(x = x, y = dspliney), size = 1.3, colour = "darkgrey") +
  geom_line(aes(x = x, y = y), size = .8) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=23),
        plot.title=element_text(hjust = .1, margin = margin(t = -4, b=-20)))




as_ggplot(gridExtra::grid.arrange(grobs = 
                                    lapply(list(dens + ggtitle("(a)"), 
                                                kern_dens + ggtitle("(b)"),
                                                spline_dens + ggtitle("(c)"), 
                                                ind_dens + ggtitle("(d)"),
                                                ord_dens + ggtitle("(e)"), 
                                                cdf_dens + ggtitle("(f)")), "+", margin), 
                                  left = text_grob("f(x)", size = 23, 
                                                   vjust = 1, rot = 90),
                                  bottom = text_grob("x", size = 23,
                                                     vjust = -.5)))






























# 
# dens_bounds %>% 
#   # pivot_longer(1:6)
#   ggplot() +
#   geom_ribbon(aes(x = x, ymin = `0.025`, ymax = `0.975`),
#               fill = "pink") +
#   geom_histogram(data = data.frame(x = samp), 
#                  aes(x = x, y = ..density..),
#                  fill = "grey") 
#   # geom_line(aes(x = x, y = `0.025`),
#   #           colour = "pink", size = 1) + 
#   # geom_line(aes(x = x, y = `0.975`),
#   #           colour = "pink", size = 1) +
#   # geom_line(aes(x = x, y = `0.25`),
#   #           colour = "red")
#   geom_ribbon(aes(x = x, ymin = `0.025`, ymax = `0.975`),
#               fill = "pink") + 
#   geom_ribbon(aes(x = x, ymin = `0.25`, ymax = `0.75`),
#               fill = "red") +
#   geom_line(aes(x = x, y = y), size = 1,
#             colour = "purple") +
#   theme_bw()



# ord_data[[2]] %>% 
#   # filter(between(prob, .05, .95)) %>%
#   ggplot() +
#   geom_segment(aes(x = prob, y = `0.025`, yend = `0.975`),
#                size = 1.5) +
#   geom_segment(aes(x = prob, y = `0.25`, yend = `0.75`), 
#                colour = "red",
#                size = 10) +
#   geom_point(aes(x = prob, y = `0.5`), colour = "pink",
#              fill = "pink", shape = 24, size = 2) +
#   geom_point(aes(x = prob, y = quantile), size = 2,
#              colour = "violet") +
#   coord_cartesian(ylim=c(-5, 5)) +
#   theme_bw()













