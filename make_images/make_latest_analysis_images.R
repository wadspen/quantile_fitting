setwd(here::here())
setwd("./make_images")
library(ggplot2)
library(dplyr)
library(stringr)
library(lubridate)


methods <- c("ind", "ord", "qgp")
scores_loc <- "../flu_forecast_analysis/"

models <- c("QGP", "ORD", "IND", "SPL", "KDE")
ltypes <- c("solid", "dotdash", "dashed", "longdash", "dotted")
colours <- c("#0072B2", "#E69F00", "#009E73",
             "#CC79A7", "#F0E442")

all_scores_og <- data.frame()
all_targets_og <- list()
for (i in 1:length(methods)) {
  scores <- readRDS(paste0(scores_loc, "all_scores_no_cover_", methods[i], ".rds")) %>% 
    mutate(method = methods[i],
           target = paste(model, date, location, horizon, sep = "_"),
           date = date(date))
  
  print(dim(scores))
  all_targets_og[[i]] <- scores$target
  all_scores_og <- rbind(all_scores_og, scores)
}

all_og_targs <- all_targets

all_scores <- readRDS("../model-fits/comb_res.rds") %>% 
  filter(str_detect(fit_model, "_sb")) %>% 
  filter(comps == 20) %>% 
  mutate(model = forc_model, method = fit_model, 
         date = date(date),
         target = paste(model, date, location, horizon, sep = "_"))



all_targets <- list()
methods <- unique(all_scores$method)
for (m in 1:length(methods)) {
   ms <- all_scores %>% 
    filter(method == methods[m])
   all_targets[[m]] <- unique(ms$target)
}

all_targs <- Reduce(intersect, all_targets)


models <- c("QGP", "ORD", "IND", "SPL", "KDE")
ltypes <- c("solid", "dotdash", "dashed", "longdash", "dotted")
colours <- c("#0072B2", "#E69F00", "#009E73",
             "#CC79A7", "#F0E442")
all_scores %>% 
  mutate(method = fit_model, 
         method = ifelse(method == "clt_sb", "QGP",
                         ifelse(method == "ord_sb", "ORD", "IND"))) %>% 
  mutate(method = factor(method, levels = c("QGP", "ORD", "IND"))) %>% 
  filter(target %in% all_targs) %>% group_by(method) %>% summarise(m = mean(fittime))
  ggplot() +
  geom_boxplot(aes(x = method, y = fittime, colour = method),
               size = 1) + 
  scale_colour_manual(name = "Model",
                      labels = models,
                      values = colours) +
  ylab("Minutes") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    panel.grid = element_blank(),
    legend.position = "none"
  )

spl_kern <- readRDS(paste0(scores_loc, "spline_kern_scores.rds"))

spl_kern <- spl_kern %>% 
  mutate(target = paste(model, date, location, horizon, sep = "_")) %>%
  mutate(date = date(date)) %>% 
  filter(target %in% targs) %>% 
  filter(zeros == TRUE) %>% 
  mutate(wis = true_wis, mae = saq, mse = ssq) %>% 
  dplyr::select(model, method, date, location, horizon, crps, logs,
                wis, est_wis, mae, mse, target)
                # true_wis, est_wis, ssq, saq, target, n, dist_diff)

all_targets_np <- list()
methods <- unique(spl_kern$method)
for (m in 1:length(methods)) {
  ms <- spl_kern %>% 
    filter(method == methods[m])
  all_targets_np[[m]] <- unique(ms$target)
}


targs_all_np <- Reduce(intersect, c(all_targets, all_targets_np))
targs_all_og <- Reduce(intersect, c(all_targets, all_targets_og))

all_scores <- all_scores %>%
  ungroup() %>% 
  filter(target %in% all_targs) %>%
  dplyr::select(model, method, date, location, horizon, crps, logs,
                wis, est_wis, mae, mse, target) 
                # true_wis0, est_wis0, ssq0, saq0, target, n0, dist_diff0)

all_scores_og <- all_scores_og %>% 
  filter(target %in% all_targets_og) %>% 
  mutate(wis = true_wis0, est_wis = est_wis0,
         mae = saq, mse = ssq) %>% 
  dplyr::select(model, method, date, location, horizon, crps, logs,
                wis, est_wis, mae, mse, target)

colnames(all_scores) <- str_replace(colnames(all_scores), "0", "")

all_scores <- rbind(all_scores, spl_kern, all_scores_og)

all_scores <- all_scores %>% 
  mutate(prior = ifelse(str_detect(method, "sb"), "DPM", "FM")) %>% 
  mutate(method = ifelse(str_detect(method, "ind"), "IND", 
                         ifelse(str_detect(method, "clt"), "QGP",
                                ifelse(str_detect(method, "ord"), "ORD",
                                       ifelse(method == "spline", "SPL",
                                              ifelse(method == "kern", "KDE",
                                                     method)))))) %>% 
  mutate(method = factor(method, levels = models))

all_cors <- all_scores %>% 
  filter(method %in% c("SPL", "KDE") | prior == "DPM") %>% 
  filter(target %in% targs_all_np) %>% 
  # filter(logs < Inf) %>% 
  filter(horizon != -1) %>% 
  mutate(horizon = paste(horizon + 1, 
                         "weeks ahead")) %>%
  group_by(method) %>% 
  summarise(
    cor = round(cor(crps, wis),3),
    corl = cor(logs, crps)
  ) 

all_scores %>% 
  filter(method %in% c("SPL", "KDE") | prior == "DPM") %>% 
  filter(target %in% targs_all_np) %>% 
  filter(horizon != -1) %>% 
  mutate(horizon = paste(horizon + 1, "weeks ahead")) %>% 
  ggplot() +
  geom_point(aes(y = crps, x = wis), alpha = .1) +
  geom_text(data = all_cors,
            aes(x = 4, y = .6, 
                label = paste0("CORR: ", cor)),
            size = 6) +
  geom_abline(slope = 1, intercept = 0, 
              colour = "darkgrey", size = .9) +
  ylab("CRPS") +
  xlab("WIS") +
  facet_wrap(~method) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17))


maep <- all_scores %>% 
  filter(method %in% c("SPL", "KDE") | prior == "DPM") %>% 
  filter(target %in% targs_all_np) %>% 
  group_by(method, date) %>% 
  summarise(msaq = mean(mae),
            mssq = mean(mse)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = msaq, colour = method, linetype = method),
            linewidth = 1.2) +
  scale_colour_manual(name = "Model",
                      labels = models,
                      values = colours) +
  scale_linetype_manual(name = "Model",
                      labels = models,
                      values = ltypes) +
  ylab("MAE") +
  xlab("") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17),
        legend.position = "none")


msep <- all_scores %>% 
  filter(method %in% c("SPL", "KDE") | prior == "DPM") %>% 
  filter(target %in% targs_all_np) %>% 
  group_by(method, date) %>% 
  summarise(msaq = mean(mae),
            mssq = mean(mse)) %>% 
  ggplot() +
  geom_line(aes(x = date, y = mssq, colour = method, linetype = method),
            linewidth = 1.2) +
  scale_colour_manual(name = "Model",
                      labels = models,
                      values = colours) +
  scale_linetype_manual(name = "Model",
                        labels = models,
                        values = ltypes) +
  ylab("MSE") +
  xlab("") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17),
        legend.position = c(.85,.8))

library(cowplot)
plot <- cowplot::plot_grid(maep, msep, nrow = 1)
ggdraw(plot) +
  draw_label("Date",
             x = 0.5, y = 0.05, size = 22,
             vjust = 0.5)

all_scores %>% 
  filter(method %in% c("QGP", "ORD", "IND")) %>% 
  ggplot() +
  geom_histogram(aes(x = method, y = fittime))



rank_scores <- all_scores %>% 
  filter(horizon != -1) %>% 
  group_by(method, model) %>% 
  summarise(
    corwc = cor(crps, true_wis),
    mcrps = mean(crps),
    mwis = mean(true_wis)
  ) %>% 
  mutate(rankc = order(mcrps),
         rankw = order(mwis)) %>% 
  mutate(ranks = rankc == rankw, rankdist = abs(rankc - rankw)) %>% 
  arrange(rankw)

rank_scores %>% 
  group_by(method) %>% 
  summarise(s = sum(ranks))


rank_scores %>% 
  ggplot() +
  geom_tile(aes(x = rankc, y = rankw, fill = corwc)) +
  facet_wrap(~method)

















