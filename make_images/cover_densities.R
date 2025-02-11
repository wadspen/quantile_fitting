source("../simulation/simulation_functions.R")
library(cmdstanr)
library(distfromq)
library(evmix)
library(dplyr)
library(tidyr)
library(ggplot2)





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




samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)

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
p <- 6

all_draws <- data.frame()
for (i in c(3, 7, 10)) {
  p <- i
  for (j in 1:4) {
    
    N <- samp_sizes[j]
    samp <- rnorm(N, mu, sigma)
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = N)
    
    
    clt_fit <- cltmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    clt_draws <- clt_fit$draws(variables = c("mu", "sigma", "n"), format = "df")
    clt_draws$model <- "clt"
    
    
    ord_fit <- ordmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    ord_draws <- ord_fit$draws(variables = c("mu", "sigma", "n"), format = "df")
    ord_draws$model <- "ord"
    
    
    cltn_fit <- cltnmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    cltn_draws <- cltn_fit$draws(variables = c("mu", "sigma"), format = "df")
    cltn_draws$n <- NA
    cltn_draws$model <- "cltn"
    
    
    ordn_fit <- ordnmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    ordn_draws <- ordn_fit$draws(variables = c("mu", "sigma"), format = "df")
    ordn_draws$n <- NA
    ordn_draws$model <- "ordn"
    
    
    ind_fit <- indmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    ind_draws <- ind_fit$draws(variables = c("mu", "sigma"), format = "df")
    ind_draws$n <- NA
    ind_draws$model <- "ind"
    
    
    
    draws <- rbind(clt_draws, ord_draws, cltn_draws, ordn_draws, ind_draws)
    draws <- draws %>% 
      mutate(N = N, quant = length(probs))
    
    all_draws <- rbind(all_draws, draws)
    
    
  }
}
# set.seed(21)
thin <- seq(1, nrow(all_draws), by = 20)
thin_draws <- all_draws[thin,]

########################
######Fit type key######
########################

models <- c("QGP", "ORD", "IND", "QGP-n", "ORD-n", "SPL", "KDE")
ltypes <- c("solid", "longdash", "dashed")
colours <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00",
             "#CC79A7", "#F0E442")

# colours <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggborderline)
mu_plot <- all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(xintercept = mu, size = .6) +
  # geom_density(aes(x = mu, group = model,
  #                  colour = model,
  #                  linetype = model), 
  #              trim = TRUE, show_guide = FALSE,
  #              size = 1.2
  #              # , linetype = c("solid", "dashed", "dotdash")
  #              ) +
  stat_density(aes(x = mu, colour = model, linetype = model),
               geom="line",position="identity", size = .7) +
  
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  # scale_colour_hue(name = "Model",
  #                  labels = c("QGP", "IND", "ORD")) +
  scale_colour_manual(name = "",
                    labels = c("QGP", "IND", "ORD")
                    ,values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_linetype_manual(name= "",
                        values=c("solid", "dashed", "longdash"),
                        labels=c("QGP", "IND","ORD")) +
  scale_x_continuous(n.breaks=3) +
  ylab("") +
  xlab(expression(mu)) +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.position = c(1,1),
    axis.text=element_text(size=8),
    axis.text.y = element_blank(),
    axis.title=element_text(size=17),
    # strip.text = element_text(size = 8),
    strip.text = element_text(margin = margin(0,0,0,0, "cm")),
    axis.ticks.y = element_blank()
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 11)
  )#, legend.position="none")

sigma_plot <- all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(xintercept = sigma, size = .6) +
  # geom_density(aes(x = mu, group = model,
  #                  colour = model,
  #                  linetype = model), 
  #              trim = TRUE, show_guide = FALSE,
  #              size = 1.2
  #              # , linetype = c("solid", "dashed", "dotdash")
  #              ) +
  stat_density(aes(x = sigma, colour = model, linetype = model),
               geom="line",position="identity", size = .7) +
  
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  # scale_colour_hue(name = "Model",
  #                  labels = c("QGP", "IND", "ORD")) +
  scale_colour_manual(name = "Model",
                      labels = c("QGP", "IND", "ORD")
                      ,values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_linetype_manual(name= "Model",
                        values=c("solid", "dashed", "longdash"),
                        labels=c("QGP", "IND","ORD")) +
  scale_x_continuous(n.breaks=3) +
  ylab("") +
  xlab(expression(sigma)) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.position = c(.1,-.4),
    axis.text=element_text(size=8),
    axis.text.y = element_blank(),
    axis.title=element_text(size=17),
    # strip.text = element_text(size = 8),
    strip.text = element_text(margin = margin(0,0,0,0, "cm")),
    axis.ticks.y = element_blank()
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 11)
  )




n_plot <- all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(xintercept = N, size = .6) +
  # geom_density(aes(x = mu, group = model,
  #                  colour = model,
  #                  linetype = model), 
  #              trim = TRUE, show_guide = FALSE,
  #              size = 1.2
  #              # , linetype = c("solid", "dashed", "dotdash")
  #              ) +
  stat_density(aes(x = n, colour = model, linetype = model),
               geom="line",position="identity", size = .7) +
  
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  # scale_colour_hue(name = "Model",
  #                  labels = c("QGP", "IND", "ORD")) +
  scale_colour_manual(name = "",
                      labels = c("QGP", "ORD")
                      ,values = c("#0072B2", "#E69F00")) +
  scale_linetype_manual(name= "",
                        values=c("solid", "longdash"),
                        labels=c("QGP","ORD")) +
  scale_x_continuous(n.breaks=3) +
  ylab("") +
  xlab("n") +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.position = c(1,1),
    axis.text=element_text(size=8),
    axis.text.y = element_blank(),
    axis.title=element_text(size=17),
    # strip.text = element_text(size = 8),
    strip.text = element_text(margin = margin(0,0,0,0, "cm")),
    axis.ticks.y = element_blank()
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 11)
    )



cowplot::plot_grid(mu_plot, sigma_plot, n_plot, nrow = 2)


all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(xintercept = sigma, size = .8) +
  # geom_density(aes(x = sigma, colour = model), trim = TRUE, show_guide = FALSE,
  #              size = 1) +
  stat_density(aes(x = sigma, colour = model, linetype = model),
               geom="line",position="identity", size = 1.2) +
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  scale_colour_hue(name = "Model",
                   labels = c("QGP", "IND", "ORD")) +
  scale_linetype_manual(name= "Model",
                        values=c("solid", "dotdash", "dashed"),
                        labels=c("QGP", "IND","ORD")) +
  ylab("") +
  xlab(expression(sigma)) +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))#, legend.position="none")




all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_density(aes(x = n, colour = model), trim = TRUE, show_guide = FALSE,
               size = 1) +
  stat_density(aes(x = n, colour = model),
               geom="line",position="identity", size = .9) +
  # geom_vline(xintercept = N, size = .9) + 
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  geom_vline(xintercept = N, size = .9) + 
  labs(colour = "Model") +
  scale_colour_hue( 
    labels = c("QGP", "ORD")) +
  ylab("") +
  xlab("N") +
  theme_bw() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        strip.text = element_text(size = 10))#, legend.position="none")


all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(aes(xintercept = N), size = .8) +
  # geom_density(aes(x = lambda, colour = model), trim = TRUE, show_guide = FALSE,
  #              size = 1) +
  stat_density(aes(x = n, colour = model, linetype = model),
               geom="line",position="identity", size = 1.2) +
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  scale_colour_hue(name = "Model",
                   labels = c("QGP", "ORD")) +
  scale_linetype_manual(name= "Model",
                        values=c("solid", "dashed"),
                        labels=c("QGP","ORD")) +
  ylab("") +
  xlab("n") +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.text.x = element_text(size=13, angle = 30),
        axis.title=element_text(size=23),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))

#############################################
################Exponential##################
#############################################


cltmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_exponential_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          '../stan_models/cdf_exponential_ind_quantiles.stan')

ordmod <- cmdstan_model(stan_file = 
                          '../stan_models/order_exponential_quantiles.stan')

cltnmod <- cmdstan_model(stan_file = 
                           '../stan_models/cdf_exponential_n_quantiles.stan')

ordnmod <- cmdstan_model(stan_file = 
                           '../stan_models/order_exponential_n_quantiles.stan')
mod_loc <- "../stan_models/"




samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)

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


lambda <- .4

all_draws <- data.frame()
for (i in c(3, 7, 10)) {
  p <- i
  for (j in 1:4) {
    
    N <- samp_sizes[j]
    samp <- rexp(N, lambda)
    probs <- levels[[p]]
    quantiles <- quantile(samp, probs)
    
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = N)
    
    
    clt_fit <- cltmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    clt_draws <- clt_fit$draws(variables = c("lambda", "n"), format = "df")
    clt_draws$model <- "clt"
    
    
    ord_fit <- ordmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    ord_draws <- ord_fit$draws(variables = c("lambda", "n"), format = "df")
    ord_draws$model <- "ord"
    
    
    cltn_fit <- cltnmod$sample(data = stan_data, 
                               iter_warmup = 10000,
                               iter_sampling = 50000,
                               chains = 1)
    
    cltn_draws <- cltn_fit$draws(variables = c("lambda"), format = "df")
    cltn_draws$n <- NA
    cltn_draws$model <- "cltn"
    
    
    ordn_fit <- ordnmod$sample(data = stan_data, 
                               iter_warmup = 10000,
                               iter_sampling = 50000,
                               chains = 1)
    
    ordn_draws <- ordn_fit$draws(variables = c("lambda"), format = "df")
    ordn_draws$n <- NA
    ordn_draws$model <- "ordn"
    
    
    ind_fit <- indmod$sample(data = stan_data, 
                             iter_warmup = 10000,
                             iter_sampling = 50000,
                             chains = 1)
    
    ind_draws <- ind_fit$draws(variables = c("lambda"), format = "df")
    ind_draws$n <- NA
    ind_draws$model <- "ind"
    
    
    
    draws <- rbind(clt_draws, ord_draws, cltn_draws, ordn_draws, ind_draws)
    draws <- draws %>% 
      mutate(N = N, quant = length(probs))
    
    all_draws <- rbind(all_draws, draws)
    
    
  }
}


lambda_plot <- all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(xintercept = lambda, size = .6) +
  # geom_density(aes(x = lambda, colour = model), trim = TRUE, show_guide = FALSE,
  #              size = 1) +
  stat_density(aes(x = lambda, colour = model, linetype = model),
               geom="line",position="identity", size = .7) +
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  scale_colour_manual(name = "Model",
                      labels = c("QGP", "IND", "ORD")
                      ,values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_linetype_manual(name= "Model",
                        values=c("solid", "dashed", "longdash"),
                        labels=c("QGP", "IND","ORD")) +
  scale_x_continuous(n.breaks=3) +
  ylab("") +
  xlab(expression(lambda)) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.position = c(.96,.95),
    axis.text=element_text(size=8),
    axis.text.y = element_blank(),
    axis.title=element_text(size=17),
    legend.key.height = unit(.3, 'cm'), 
    legend.key.width = unit(.4, 'cm'),
    # strip.text = element_text(size = 8),
    strip.text = element_text(margin = margin(0,0,0,0, "cm")),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7)
  )#, legend.position="none")




n_plot <- all_draws %>% 
  # filter(quant == 15, n == 50) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  filter(quant %in% c(7, 15, 23)) %>% 
  ggplot() + 
  geom_vline(aes(xintercept = N), size = .6) +
  # geom_density(aes(x = lambda, colour = model), trim = TRUE, show_guide = FALSE,
  #              size = 1) +
  stat_density(aes(x = n, colour = model, linetype = model),
               geom="line",position="identity", size = .7) +
  
  # coord_cartesian(xlim = c(3.3, 3.8)) +
  ggh4x::facet_grid2(N~quant, 
                     scales = "free", independent = "all") +
  # labs(colour = "Model") +
  scale_colour_manual(name = "",
                      labels = c("QGP", "ORD")
                      ,values = c("#0072B2", "#E69F00")) +
  scale_linetype_manual(name= "",
                        values=c("solid", "longdash"),
                        labels=c("QGP","ORD")) +
  scale_x_continuous(n.breaks=3) +
  ylab("") +
  xlab("n") +
  # ggtitle("yourmom") +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.position = c(1,1),
    axis.text=element_text(size=8),
    axis.text.y = element_blank(),
    axis.title=element_text(size=17),
    # strip.text = element_text(size = 8),
    strip.text = element_text(margin = margin(0,0,0,0, "cm")),
    axis.ticks.y = element_blank()
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 11)
  )#, legend.position="none")



cowplot::plot_grid(lambda_plot, n_plot, nrow = 1)


                      