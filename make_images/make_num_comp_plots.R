library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(kableExtra)

score_loc <- "../simulation/sim_scores/"
file_loc <- "_numcomp/size1000_probs23_scores.csv"

dists <- c("lp", "evd", "gmix")
methods <- c("MCMC", "variational")
samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
probs_sizes <- 23
all_comps <- data.frame()

for (i in 1:length(dists)) {
  for (j in 1:length(methods)) {
    for (k in 1:length(samp_sizes)) {
      for (l in 1:length(probs_sizes)) {
        file_path <- paste0(score_loc, dists[i], "_", methods[j], 
                            "_numcomp/size", samp_sizes[k], "_probs",
                            probs_sizes[l], "_scores.csv")
        
        comps <- read.csv(file_path)
        comps$dist <- dists[i]
        comps$method <- methods[j]
        all_comps <- rbind(all_comps, comps)
        
      }
    }
  }
}


sum_scores <- all_comps %>% 
  group_by(n, probs, quants, model, method, dist) %>% 
  summarise(muwd1 = mean(uwd1, na.rm = TRUE),
            muwd2 = mean(uwd2, na.rm = TRUE),
            mks = mean(ks, na.rm = TRUE),
            mkld = mean(kld, na.rm = TRUE),
            mtv = mean(tv, na.rm = TRUE)) %>% 
  mutate(comps = as.numeric(str_extract(model, "\\d+"))) %>% 
  mutate(n_known = str_detect(model, "n")) %>% 
  ungroup() %>% 
  select(-probs, -quants)


uwdp <- sum_scores %>% 
  filter(method == "MCMC", n == 1000) %>%  
  mutate(dist = ifelse(dist == "evd", "EV",
                       ifelse(dist == "lp", "La", 
                              ifelse(dist == "gmix", "MIX", dist)))) %>% 
  mutate(n_known = ifelse(n_known == TRUE, "QGP-n", "QGP")) %>% 
  ggplot() +
  geom_line(aes(x = comps, y = muwd1, colour = dist,
                linetype = dist), size = 1.5) +
  facet_grid(~n_known, scales = "free") +
  scale_linetype_manual(values = c("EV" = "solid", "La" = "longdash", 
                                   "MIX" = "dashed")) +
  scale_colour_manual(values = c("EV" = "#0072B2", "La" = "#E69F00", 
                                 "MIX" = "#009E73")) +
  scale_x_continuous(breaks = 1:6) +
  xlab("Number of Mixture Distribution Components") +
  ylab("UWD1") +
  theme_bw() +
  theme(legend.position = c(.88, .65),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.title.x =element_blank(),
        strip.text = element_text(size = 12))



tvp <- sum_scores %>% 
  filter(method == "variational", n == 1000) %>%  
  mutate(dist = ifelse(dist == "evd", "EV",
                       ifelse(dist == "lp", "La", 
                              ifelse(dist == "gmix", "MIX", dist)))) %>% 
  mutate(n_known = ifelse(n_known == TRUE, "QGP-n", "QGP")) %>% 
  ggplot() +
  geom_line(aes(x = comps, y = mtv, colour = dist,
                linetype = dist), size = 1.5) +
  facet_grid(~n_known, scales = "free") +
  scale_linetype_manual(values = c("EV" = "solid", "La" = "longdash", 
                                   "MIX" = "dashed")) +
  scale_colour_manual(values = c("EV" = "#0072B2", "La" = "#E69F00", 
                                 "MIX" = "#009E73")) +
  scale_x_continuous(breaks = 1:6) +
  xlab("Number of Mixture Distribution Components") +
  ylab("TV") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(l = 30, ),
        strip.text = element_blank())

kldp <- sum_scores %>% 
  filter(method == "variational", n == 1000) %>%  
  mutate(dist = ifelse(dist == "evd", "EV",
                       ifelse(dist == "lp", "La", 
                              ifelse(dist == "gmix", "MIX", dist)))) %>% 
  mutate(n_known = ifelse(n_known == TRUE, "QGP-n", "QGP")) %>% 
  ggplot() +
  geom_line(aes(x = comps, y = mkld, colour = dist,
                linetype = dist), size = 1.5) +
  facet_grid(~n_known, scales = "free") +
  scale_linetype_manual(values = c("EV" = "solid", "La" = "longdash", 
                                   "MIX" = "dashed")) +
  scale_colour_manual(values = c("EV" = "#0072B2", "La" = "#E69F00", 
                                 "MIX" = "#009E73")) +
  scale_x_continuous(breaks = 1:6) +
  xlab("Number of Mixture Distribution Components") +
  ylab("KLD") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 14),
        strip.text = element_blank())



library(patchwork)
uwdp / tvp / kldp


sum_scores %>%
  filter(method == "MCMC", n == 1000, n_known == TRUE, dist == "evd") %>% 
  arrange(comps) %>% 
  ungroup() %>% 
  mutate(eff = lag(muwd1)/muwd1)

uwd_sum <- sum_scores %>% 
  filter(n == 1000, method == "MCMC") %>% 
  mutate(n_known = ifelse(n_known == TRUE, "QGP-n", "QGP")) %>%
  select(n_known, dist, comps, muwd1) %>% 
  pivot_wider(names_from = comps, values_from = muwd1)

sum_eff <- sum_scores %>% 
  filter(n == 1000) %>% 
  arrange(method, n_known, dist, n, comps) %>% 
  group_by(method, n_known, dist, n) %>% 
  mutate(eff_uwd1 = lag(muwd1)/muwd1,
         eff_kld = lag(mkld)/mkld,
         eff_tv = lag(mtv)/mtv)



eff_table <- sum_eff %>% 
  filter(!is.na(eff_uwd1), method == "MCMC", n_known == TRUE, n == 1000) %>% 
  mutate(comps = as.factor(paste0("c", comps))) %>% 
  ungroup() %>% 
  select(dist, comps, eff_uwd1) %>% 
  # mutate(eff_uwd1 = (eff_uwd1 - 1) * 100) %>% 
  pivot_wider(
    names_from = comps,
    values_from = eff_uwd1
  )


eff_table %>%
  kable("latex", booktabs = TRUE, align = "l l c c c c c") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "top")


scores_sum <- sum_scores %>% 
  filter(n == 1000, method == "variational") %>% 
  mutate(n_known = ifelse(n_known == TRUE, "QGP-n", "QGP")) %>%
  select(n_known, dist, comps, mkld) %>% 
  mutate(mkld = round(mkld, 3)) %>% 
  pivot_wider(names_from = comps, values_from = mkld)


scores_sum %>% 
  kable("latex", booktabs = TRUE, align = "l l c c c c c c") %>%
  kable_styling(latex_options = c("hold_position"))
  




















