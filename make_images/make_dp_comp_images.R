setwd(here::here())
setwd("./make_images/")

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

ncomps <- c(2,3,4,6,8,12,20)

bayes_res <- data.frame()
for (i in 1:length(ncomps)) {
  res <- readRDS(paste0("../simulation/simple_res_tm/ncomp",
                        ncomps[i], "/comb_res.rds"))
  res$ncomp <- ncomps[i]
  
  bayes_res <- rbind(bayes_res, res)
}

# mods <- unique(bayes_res$model[bayes_res$ncomp == 2])
# ns <- unique(bayes_res$n[bayes_res$ncomp == 2])

freq_res <- readRDS("../simulation/simple_res/comb_res_freq.rds") %>% 
  mutate(ncomp = NA)

# bayes_res <- readRDS("../simulation/simple_res_tm/ncomp20/comb_res.rds") %>% 
#   mutate(ncomp = 20)

all_res <- rbind(bayes_res, freq_res) %>% 
  filter(dist != "tdp") %>% 
  # filter(n %in% ns, model %in% c(mods, "spline", "kern"), dist != "tdp") %>%
  mutate(sb = ifelse(str_detect(model, "_sb"), "DPM", "FM")) %>% 
  mutate(mod = str_remove(model, "_sb")) %>% 
  mutate(dist = ifelse(dist == "evd", "EV",
                       ifelse(dist == "lp", "La", "MIX"))) %>% 
  mutate(mod = ifelse(mod == "clt", "QGP", 
                      ifelse(mod == "ind", "IND", 
                             ifelse(mod == "ord", "ORD", 
                                    ifelse(mod == "spline", "SPL", 
                                           ifelse(mod == "kern", "KDE", NA)))))) %>% 
  mutate(mod = factor(mod, levels = c("QGP", "IND","KDE", "ORD", "SPL")))



res_sum <- all_res %>% 
  group_by(mod, sb, n, numq, dist, ncomp) %>% 
  summarise(`UWD1` = mean(uwd1),
            TV = mean(tv),
            KLD = mean(kld),
            mtime = median(ftime),
            `% Coverage` = mean(coverm), 
            `Est. Comp` = mean(comp99)) 

res_sum %>% 
  filter(n == 500, sb == "DPM") %>%
  ggplot() +
  geom_line(aes(x = numq, y = mtime, colour = mod, linetype = mod)) +
  facet_grid(ncomp~dist, scale = "free_y")  

all_res %>% 
  filter(sb == 'DPM', ncomp == 20, numq %in% c(11, 15, 23),
         n %in% c(500)) %>% 
  ggplot() +
  geom_boxplot(aes(x = mod, y = comp99)) +
  ylab("Estimated # Mixture Components") +
  xlab("") +
  facet_grid(dist~numq) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 19),
    strip.text = element_text(size = 17)
  )

res_sum %>% 
  filter(sb == "DPM", ncomp == 20, numq %in% c(5, 11, 19, 23)) %>% 
  ggplot() +
  geom_line(aes(x = n, y = `Est. Comp`, colour = mod, linetype = mod)) +
  facet_grid(numq ~ dist)



res_sum %>% 
  filter(numq %in% 11:23, ncomp %in% 8:20, dist == "MIX", n > 50) %>% 
  # filter(ncomp == 20, numq == 23) %>% 
  filter(mod %in% c("QGP", "ORD", "IND")) %>% 
  ggplot() +
  geom_line(aes(x = n, y = `UWD1`, colour = sb, linetype = mod), size = 1.3) +
  facet_grid(ncomp ~ numq, scales = "free_y") +
  plot_annotation(
    title = "# Quantiles K",
    theme = theme(
      plot.title = element_text(
        # angle = -90,
        vjust = 0.5,
        hjust = 0.5,
      )
    )
  ) +
  scale_colour_manual(name = "Prior",
                      labels = c("DPM","FM"),
                      values = c("#0072B2", "#E69F00")) +
  scale_linetype_manual(name= "Model",
                        values=c(1:2,4),
                        labels=c("QGP", "IND", "ORD")) +
  ylab("UWD1") +
  xlab("n") +
  scale_x_continuous(n.breaks = 3) +
  theme_bw() +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size = 12, angle = 0, hjust =1),
        axis.title=element_text(size=16),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom") +
  wrap_elements(
    full = textGrob(
      "# Normal Components",
      rot = 270,
      gp = gpar(fontsize = 16)
    )
  ) +
  plot_layout(widths = c(1, 0.06))


res_sum %>% 
  filter((sb == "FM" & ncomp == 8) | (sb == "DPM" & ncomp == 20)) %>% 
  filter(n > 50) %>%
  ggplot() +
  geom_line(aes(x = numq, y = `UWD1`, linetype = mod, colour = sb),
            size = 1.3) +
  facet_grid(dist~n, scales = "free_y") +
  scale_colour_manual(name = "Prior",
                      labels = c("DPM","FM"),
                      values = c("#0072B2", "#E69F00")) +
  scale_linetype_manual(name= "Model",
                        values=c(1:2,4),
                        labels=c("QGP", "IND", "ORD")) +
  ylab("UWD1") +
  xlab("# Quantiles K") +
  scale_x_continuous(n.breaks = 5) +
  theme_bw() +
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(size = 12, angle = 0, hjust =1),
        axis.title=element_text(size=16),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom") #+
  # wrap_elements(
  #   full = textGrob(
  #     "# Normal Components",
  #     rot = 270,
  #     gp = gpar(fontsize = 16)
  #   )
  # ) +
  # plot_layout(widths = c(1, 0.06))

  
res_sum %>% 
  filter(sb == "DPM" | mod %in% c("SPL", "KDE")) %>% 
  filter(numq == 23) %>% 
  filter(ncomp == 20 | is.na(ncomp)) %>% 
  dplyr::select(-mtime) %>% 
  pivot_longer(7:9, names_to = "meas", values_to = "mdist") %>% 
  # filter(n > 50) %>% 
  # filter(mod == "ORD") %>% 
  ggplot() +
  geom_line(aes(x = n, y = mdist, colour = mod, linetype = mod),
            size = 1.1) +
  facet_grid(meas ~ dist, scales = "free_y") +
  scale_colour_manual(name = "Model",
                      labels = c("QGP", "IND","KDE", "ORD", "SPL"),
                      values = c("#0072B2", "#009E73", "darkgrey", "#E69F00",
                                 "#CC79A7")) +
  scale_linetype_manual(name= "Model",
                        values=1:5,
                        labels=c("QGP", "IND","KDE", "ORD", "SPL")) +
  ylab("") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12, angle = 0, hjust =1),
        axis.title=element_text(size=16),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))


res_sum %>% 
  filter(sb == "DPM") %>% 
  filter(ncomp == 20, numq %in% c(3, 11, 23)) %>% 
  ggplot() +
  geom_line(aes(x = n, y = `% Coverage`*100, colour = mod, linetype = mod),
            size = 1.1) +
  geom_hline(aes(yintercept = 95), size = .8) +
  facet_grid(numq~dist, scales = "free_y") +
  scale_colour_manual(name = "Model",
                      labels = c("QGP", "IND","ORD"),
                      values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_linetype_manual(name= "Model",
                        values=c(1,2,4),
                        labels=c("QGP", "IND","ORD")) +
  ylab("% Coverage Q(p)") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12, angle = 0, hjust =1),
        axis.title=element_text(size=16),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))





res_sum %>% 
  filter(sb %in% c("FM")) %>% 
  filter(numq %in% c(11, 15, 23), mod == "QGP", n == 500) %>% 
  # dplyr::select(-mtime) %>% 
  # pivot_longer(7:9, names_to = "meas", values_to = "mdist") %>%
  ggplot() +
  geom_line(aes(x = ncomp, y = `UWD1`, colour = dist, linetype = dist),
            size = 1.1) +
  scale_colour_manual(name = "",
                      labels = c("EV", "La","MIX"),
                      values = c("#0072B2", "#009E73", "#E69F00")) +
  scale_linetype_manual(name= "",
                        values=c(1,2,4),
                        labels=c("EV", "La","MIX")) +
  # ylab("% Coverage") +
  facet_grid(numq~., scales = "free_y") +
  scale_y_continuous(n.breaks = 3) +
  xlab("# Quantiles K") +
  theme_bw() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x=element_text(size = 12, angle = 0, hjust =1),
        axis.title=element_text(size=16),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.position = c(.2,.87))




