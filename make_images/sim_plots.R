library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

all_sum_scores <- data.frame()
for (type in c("VI", "MCMC")) {
  for (true in c("gmix", "evd", "lp")) {
    name <- paste0("../simulation/comb_res/scores_", true, "_", type, ".rds")
    score <- readRDS(name)
    
   
    true <- ifelse(true == "gmix", "MIX", 
                   ifelse(true == "lp", "La", "EV"))
    sum_score <- score %>% 
      group_by(n, quants, model) %>% 
      summarise(mwd1 = mean(uwd1),
                sdwd1 = sd(uwd1),
                mtv = mean(tv),
                mkld = mean(kld)) %>% 
      mutate(upp95 = mwd1 + 1.96*sdwd1/sqrt(n), 
             low95 = mwd1 - 1.96*sdwd1/sqrt(n)) %>% 
      mutate(tru_mod = true)
    
    sum_score$samp_type <- type 
    
    all_sum_scores <- rbind(all_sum_scores, sum_score)
  }
  
}

models <- c("QGP", "ORD", "IND", "QGP-n", "ORD-n", "SPL", "KDE")
ltypes <- c("solid", "longdash", "dashed")
colours <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00",
             "#CC79A7", "#F0E442")

dist_plot <- all_sum_scores %>% 
  filter(quants %in% c(9,13,23,50)) %>%
  filter(quants == 23) %>% 
  filter(!(model %in% c("cltn", "ordn", "meta"))) %>%
  pivot_longer(cols = c("mwd1", "mtv", "mkld"), names_to = "measure",
               values_to = "dist") %>% 
  filter((measure %in% c("mtv", "mkld") & samp_type == "VI") |
           (measure == "mwd1" & samp_type == "MCMC")) %>% 
  mutate(measure = ifelse(measure == "mwd1", "UWD1",
                          ifelse(measure == "mtv", "TV", "KLD"))) %>% 
  ggplot() +
  geom_path(aes(x = n, y = dist, group = model, 
                colour = model, linetype = model), size = .8, alpha = .9) +
  facet_grid(measure~tru_mod, scale = "free") +
  scale_colour_manual(name = "Model",
                   labels = c("QGP", "IND","KDE", "ORD", "SPL"),
                   values = c("#0072B2", "#009E73", "darkgrey", "#E69F00",
                              "#CC79A7")) +
  scale_linetype_manual(name= "Model",
                        values=1:5,
                        labels=c("QGP", "IND","KDE", "ORD", "SPL")) +
  ylab("") +
  xlab("n") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12, angle = 30),
        axis.title=element_text(size=14),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11))
  
  
  
  
  
# cover <- readRDS("./simulation/comb_res/coverage_evd.rds")
all_cover <- data.frame()
for (true in c("gmix", "evd", "lp")) {
  name <- paste0("../simulation/comb_res/coverage_", true, ".rds")
  cover <- readRDS(name)
  
  
  true <- ifelse(true == "gmix", "MIX", 
                 ifelse(true == "lp", "La", "EV"))
  sum_cover <- cover %>% 
    group_by(n, quants, model) %>% 
    summarise(mc95 = mean(cover95),
              sdc95 = sd(cover95),
              mc50 = mean(cover50),
              sdc50 = sd(cover50)) %>% 
    # mutate(upp95 = mwd1 + 1.96*sdwd1/sqrt(n), 
    #        low95 = mwd1 - 1.96*sdwd1/sqrt(n)) %>% 
    mutate(tru_mod = true)
  
  all_cover <- rbind(all_cover, sum_cover)
}




cover_plot <- all_cover %>% 
  mutate(quants = as.numeric(quants), n = as.numeric(n)) %>% 
  filter(quants %in% c(9,13,23,50)) %>%
  filter(quants == 23) %>% 
    filter(!(model %in% c("cltn", "ordn", "meta"))) %>% 
    # filter(n > 500) %>%  
    ggplot() +
    geom_hline(yintercept = .95*100, size = .72) +
    geom_line(aes(x = n, y = mc95*100, group = model, 
                  colour = model, linetype = model), size = .8, alpha = .9) +
    # facet_grid(quants~tru_mod, scale = "free") +
  facet_wrap(~tru_mod) +
    # geom_hline(yintercept = .5, size = .72) +
  scale_x_continuous(breaks=seq(0, 5000, 2500)) +
  scale_colour_manual(name = "Model",
                   labels = c("QGP", "IND", "ORD")
                   ,values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_linetype_manual(name= "Model",
                        values=c("solid", "dashed", "longdash"),
                        labels=c("QGP", "IND","ORD")) +
    ylab("% Coverage") +
    xlab("") +
    theme_bw() +
    theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12, angle = 30),
        axis.title=element_text(size=14),
        strip.text.y = element_text(size = 12,),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11))


cowplot::plot_grid(dist_plot, cover_plot, nrow = 2, rel_heights = c(2,1))


all_cover %>% 
  filter(quants %in% c(9,13,23,50)) %>%
  mutate(quants = as.numeric(quants), n = as.numeric(n)) %>%
  filter(!(model %in% c("cltn", "ordn", "meta"))) %>% 
  # filter(n > 500) %>%  
  ggplot() +
  geom_hline(yintercept = .95, size = .72) +
  geom_line(aes(x = n, y = mc95, group = model, 
                colour = model, linetype = model), size = 1.3) +
  facet_grid(quants~tru_mod, scale = "free") +
  # geom_hline(yintercept = .95, size = .72) +
  scale_x_continuous(breaks=seq(0, 5000, 2500)) +
  scale_colour_hue(name = "Model",
                   labels = c("QGP", "IND", "ORD")) +
  scale_linetype_manual(name= "Model",
                        values=c("solid", "dotdash", "dashed"),
                        labels=c("QGP", "IND","ORD")) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12, angle = 40),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 16),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))





sum_cover %>% 
  # mutate(n = as.numeric(n)) %>% 
  # mutate(quants = as.numeric(quants)) %>% 
  mutate(prob = as.numeric(prob)) %>% 
  filter(model != "ind") %>% 
  ggplot() +
  geom_line(aes(x = prob, y = m95, colour = model)) +
  facet_grid(prob~model)



#########################################
###############Exp Scores################
#########################################

exp_score <- readRDS("./simulation/comb_res/scores_exp_MCMC.rds")


exp_sum <- exp_score %>% 
  group_by(model, n, quants) %>% 
  filter(!(model %in% c("spline", "kern", "ind"))) %>%
  # filter(n > 150) %>% 
  summarise(mwd1 = mean(uwd1),
            mwd2 =  mean(uwd2),
            mtv = mean(tv),
            mkld = mean(kld),
            mcl = mean(cover90_lambda),
            mcn = mean(cover90_n),
            mwl = mean(width_lambda),
            mwn = mean(width_n),
            mt = mean(time)*60) 

exp_sum %>% 
  filter(n %in% c(50, 150, 500, 1000, 5000), quants %in% c(3,7,15,23,50)) %>% 
  dplyr::select(model, n, quants, mcl, mcn) %>% 
  pivot_longer(4:5, names_to = "param", values_to = "mean_param") %>% 
  mutate(facet = ifelse(param == "mcl", "lambda", "n")) %>%
  # filter(param == "mcl") %>% 
  ggplot() +
  geom_hline(yintercept = 90, size = .9) +
  geom_path(aes(x = n, y = mean_param*100, group = model, 
                colour = model, linetype = model), size = 1.3) +
  # geom_path(aes(x = n, y = mcs, linetype = model, colour = model), 
  #           size = 1.3) +
  facet_grid(quants~facet, labeller = label_parsed, scale = "free") +
  xlab("n") +
  ylab("% Coverage") +
  # expand_limits(y=100) +
  ylim(c(65,100)) +
  # labs(colour = "Model") +
  scale_colour_hue(name = "Model", 
                   labels = c("QGP", "QGPN","IND", "ORD", "ORDN")) +
  scale_linetype_manual(name= "Model", 
                        values=c("solid", "dashed", "dotdash", 
                                 "solid", "dashed"),
                        labels = c("QGP", "QGPN","IND", "ORD", "ORDN")) +
  # guides(linetype = FALSE) +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12, angle=30),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 22),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))


exp_dist <- exp_score %>% 
  mutate(n = factor(n)) %>%
  group_by(model, n, quants) %>% 
  # filter(!(model %in% c("spline", "kern"))) %>%
  # filter(n > 500) %>% 
  summarise(mwd1 = mean(uwd1),
            mwd2 =  mean(uwd2),
            mkld = mean(kld),
            mtv = mean(tv)) %>% 
  pivot_longer(4:7, names_to = "metric", values_to = "dist")


exp_dist %>% 
  filter(quants %in% c(3,7,15,23,50), metric != "mwd2") %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  mutate(facet = ifelse(metric == "mwd1", "UWD1", 
                        ifelse(metric == "mtv", "TV", "KLD"))) %>% 
  mutate(facet = factor(facet, levels = c("UWD1", "TV", "KLD"))) %>%  
  ggplot() +
  geom_path(aes(x = n, y = dist, group = model, 
                colour = model, linetype = model), size = 1.3) +
  facet_grid(quants~facet, scale = "free") +
  scale_colour_hue(name = "Model", 
                   labels = c("QGP", "IND","KDE", "ORD", "SPL")) +
  scale_linetype_manual(name= "Model", 
                        values=1:5,
                        labels=c("QGP", "IND","KDE", "ORD", "SPL")) +
  ylab("") +
  xlab("n") +
  # facet_wrap(~metric) +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12, angle=30),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 18),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))









#########################################
###############Norm Scores###############
#########################################

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
setwd("./make_images/")
norm_score <- readRDS("../simulation/comb_res/scores_norm_MCMC.rds")


norm_sum <- norm_score %>% 
  mutate(n = factor(n)) %>%
  group_by(model, n, quants) %>% 
  filter(!(model %in% c("spline", "kern"))) %>%
  # filter(n > 500) %>% 
  summarise(mwd1 = mean(uwd1),
            mwd2 =  mean(uwd2),
            mkld = mean(kld),
            mtv = mean(tv),
            mcm = mean(cover90_mu),
            mcs = mean(cover90_sigma),
            mcn = mean(cover90_n),
            mwm = mean(width_mu),
            mws = mean(width_sigma),
            mwn = mean(width_n),
            mt = mean(time)*60) 



models <- c("QGP", "ORD", "IND", "QGP-n", "ORD-n", "SPL", "KDE")
ltypes <- c("solid", "longdash", "dashed")
colours <- c("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00",
             "#CC79A7", "#F0E442")


cover_plot <- norm_sum %>% 
  filter(n %in% c(50, 150, 500, 1000, 5000), quants %in% c(3,7,15,23,50)) %>% 
  filter(quants == 23) %>% 
  dplyr::select(model, n, quants, mcm, mcs, mcn) %>% 
  pivot_longer(4:6, names_to = "param", values_to = "mean_param") %>% 
  # filter(param != "mcn") %>% 
  mutate(facet = ifelse(param == "mcm", "mu", 
                        ifelse(param == "mcs", "sigma", "n"))) %>% 
  mutate(facet = factor(facet, levels = c("mu", "sigma", "n"))) %>% 
  ggplot() +
  geom_hline(yintercept = 90, size = .7) +
  geom_path(aes(x = n, y = mean_param*100, group = model, 
                colour = model, linetype = model), size = .7, alpha = .9, 
            lineend = "round") +
  # scale_x_discrete(labels=c('label1', 'label2', 'label3', 'label4', 'label5')) +
  # geom_path(aes(x = n, y = mcs, linetype = model, colour = model), 
  #           size = 1.3) +
  facet_grid(quants~facet, labeller = label_parsed, scale = "free") +
  xlab("n") +
  ylab("% Coverage") +
  # xlim(c(50, 5000)) +
  # expand_limits(y=100) +
  ylim(c(10,100)) +
  # labs(colour = "Model") +
  scale_colour_manual(name = "Model", 
                   labels = c("QGP", "QGP-n","IND", "ORD", "ORD-n"),
                   values = c("#0072B2", "#56B4E9", "#009E73",
                              "#E69F00", "#D55E00")) +
  scale_linetype_manual(name= "Model", 
                        values=c("solid", "dashed", "dotdash", 
                                 "solid", "dashed"),
                        labels = c("QGP", "QGP-n","IND", "ORD", "ORD-n")) +
  # guides(linetype = FALSE) +
  theme_bw() +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10, angle=35),
        axis.title=element_text(size=15),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 15),
        legend.title = element_text(size = 13),
        strip.text = element_text(margin = margin(0,0,0,0, "cm")),
        legend.text = element_text(size = 11),
        
        legend.key.height = unit(.3, 'cm'), 
        legend.key.width = unit(.4, 'cm'))




norm_dist <- norm_score %>% 
  filter(kld < Inf) %>% 
  mutate(n = factor(n)) %>%
  group_by(model, n, quants) %>% 
  # filter(!(model %in% c("spline", "kern"))) %>%
  # filter(n > 500) %>% 
  summarise(mwd1 = mean(uwd1),
            mwd2 =  mean(uwd2),
            mkld = mean(kld),
            mtv = mean(tv)) %>% 
  pivot_longer(4:7, names_to = "metric", values_to = "dist")


dist_plot <- norm_dist %>% 
  filter(quants %in% c(3,7,15,23,50), metric != "mwd2") %>% 
  filter(quants == 23) %>% 
  filter(!(model %in% c("cltn", "ordn"))) %>% 
  mutate(facet = ifelse(metric == "mwd1", "UWD1", 
                        ifelse(metric == "mtv", "TV", "KLD"))) %>% 
  mutate(facet = factor(facet, levels = c("UWD1", "TV", "KLD"))) %>% 
  ggplot() +
  geom_path(aes(x = n, y = dist, group = model, 
                colour = model, linetype = model), size = .7, alpha = .9) +
  facet_grid(quants~facet, scale = "free") +
  scale_colour_manual(name = "Model", 
                   labels = c("QGP", "IND","KDE", "ORD", "SPL"),
                   values = c("#0072B2", "#009E73", "darkgrey",
                     "#E69F00", "#CC79A7")) +
  scale_linetype_manual(name= "Model",
                        labels=c("QGP", "IND","KDE", "ORD", "SPL"),
                        values=c("solid", "dashed", "dotdash", 
                                 "longdash", "dotdash")) +
  scale_y_continuous(n.breaks = 3) +
  ylab("") +
  xlab("n") +
  # facet_wrap(~metric) +
  theme_bw() +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10, angle=35),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        legend.title = element_text(size = 13),
        strip.text = element_text(margin = margin(0,0,0,0, "cm")),
        legend.text = element_text(size = 11),
        
        legend.key.height = unit(.3, 'cm'), 
        legend.key.width = unit(.4, 'cm'))

cowplot::plot_grid(cover_plot, dist_plot, ncol = 1,
                   align = "hv")
  

norm_score %>% 
  group_by(model, n, quants) %>% 
  # filter(!(model %in% c("spline", "kern"))) %>% 
  filter(!(model %in% c("clt", "ordn"))) %>%
  summarise(mt = mean(time)*60) %>% 
  ggplot() +
  geom_line(aes(x = n, y = mt, colour = model)) +
  facet_wrap(~quants, scale = "free")
  



#############################################
############Tukey Lambda Time################
#############################################


##################################################
######Need to fix quant values for Tukey sim######
##################################################

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
quants <- unlist(lapply(levels, FUN = length))
fixdf <- data.frame(probs = 1:length(quants), quant = quants)

tuk_score <- readRDS("../simulation/comb_res/scores_tuk.rds")


tuk_score <- tuk_score %>% 
  left_join(fixdf, by = "probs")


tuk_score %>% 
  filter(quant %in% c(3,7,15,23,50)) %>%
  # group_by(model, quant, n) %>% 
  # summarise(mt = mean(time)*60) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(quant), y = time*60, colour = model)) +
  # scale_y_continuous(trans = "log10")
  scale_y_log10() +
  scale_colour_hue(name = "Model", 
                   labels = c("QGP", "ORD")) +
  ylab("Seconds") +
  xlab("K") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=18),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 22),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))
  
  geom_line(aes(x = n, y = mt, colour = model)) +
  facet_wrap(~quant, scale = "free")


tuk_score %>% 
  filter(probs %in% c(3,7,15,23)) %>% 
  group_by(model, n, probs) %>% 
  summarise(mcl = mean(cover90_lambda),
            mcn = mean(cover90_n),
            mwl = mean(width_lambda),
            mwn = mean(width_n)) %>% 
  ggplot() +
  geom_line(aes(x = n, y = mcn, colour = model)) +
  facet_wrap(~probs, scale = "free")
  


