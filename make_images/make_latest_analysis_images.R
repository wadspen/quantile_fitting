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

all_scores <- data.frame()
all_targets <- list()
for (i in 1:length(methods)) {
  scores <- readRDS(paste0(scores_loc, "all_scores_no_cover_", methods[i], ".rds")) %>% 
    mutate(method = methods[i],
           target = paste(model, date, location, horizon, sep = "_"),
           date = date(date))
  
  print(dim(scores))
  all_targets[[i]] <- scores$target
  all_scores <- rbind(all_scores, scores)
}

targs <- Reduce(intersect, all_targets)

spl_kern <- readRDS(paste0(scores_loc, "spline_kern_scores.rds"))

spl_kern <- spl_kern %>% 
  mutate(target = paste(model, date, location, horizon, sep = "_")) %>%
  mutate(date = date(date)) %>% 
  filter(target %in% targs) %>% 
  filter(zeros == TRUE) %>% 
  dplyr::select(model, method, date, location, horizon, crps, logs,
                true_wis, est_wis, ssq, saq, target, n, dist_diff)

all_scores <- all_scores %>%
  filter(target %in% targs) %>%
  dplyr::select(model, method, date, location, horizon, crps, logs,
                true_wis0, est_wis0, ssq0, saq0, target, n0, dist_diff0)

colnames(all_scores) <- str_replace(colnames(all_scores), "0", "")

all_scores <- rbind(all_scores, spl_kern)

all_scores <- all_scores %>% 
  mutate(method = ifelse(method == "ind", "IND", 
                         ifelse(method == "qgp", "QGP",
                                ifelse(method == "ord", "ORD",
                                       ifelse(method == "spline", "SPL",
                                              ifelse(method == "kern", "KDE",
                                                     method)))))) %>% 
  mutate(method = factor(method, levels = models))

all_cors <- all_scores %>% 
  # filter(logs < Inf) %>% 
  filter(horizon != -1) %>% 
  mutate(horizon = paste(horizon + 1, 
                         "weeks ahead")) %>%
  group_by(method) %>% 
  summarise(
    cor = round(cor(crps, true_wis),3),
    corl = cor(logs, crps)
  ) 

all_scores %>% 
  filter(horizon != -1) %>% 
  mutate(horizon = paste(horizon + 1, "weeks ahead")) %>% 
  ggplot() +
  geom_point(aes(y = crps, x = true_wis), alpha = .1) +
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
  group_by(method, date) %>% 
  summarise(msaq = mean(saq),
            mssq = mean(ssq)) %>% 
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
  xlab("Date") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17),
        legend.position = "none")


msep <- all_scores %>% 
  group_by(method, date) %>% 
  summarise(msaq = mean(saq),
            mssq = mean(ssq)) %>% 
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
  xlab("Date") +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17),
        legend.position = c(.85,.8))


cowplot::plot_grid(maep, msep, nrow = 1)

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

















