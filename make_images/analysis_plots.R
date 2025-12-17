library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)



all_scores <- readRDS("../flu_forecast_analysis/all_scores.rds")

all_cors <- all_scores %>% 
  # filter(logs < Inf) %>% 
  filter(horizon != -1) %>% 
  mutate(horizon = paste(horizon + 1, 
                         "weeks ahead")) %>%
  group_by(horizon) %>% 
  summarise(
    cor = round(cor(crps, true_wis),2),
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
  facet_wrap(~horizon) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17))


all_scores %>% 
  group_by(model, horizon) %>% 
  summarise(
    mcrps = mean(crps),
    mwis = mean(true_wis)
  ) %>% 
  ggplot() +
  geom_point(aes(x = mcrps, y = mwis)) +
  facet_wrap(~horizon)

all_scores %>% 
  filter(horizon != -1) %>% 
  mutate(date = date(date)) %>% 
  group_by(date, horizon) %>% 
  summarise(
    cor = cor(crps, true_wis)
  ) %>% 
  ggplot() +
  geom_point(aes(x = date, y = cor)) +
  facet_wrap(~horizon)
  

rank_scores <- all_scores %>% 
  filter(horizon != -1) %>% 
  group_by(model) %>% 
  summarise(
    corwc = cor(crps, true_wis),
    mcrps = mean(crps),
    mwis = mean(true_wis)
  ) %>% 
  mutate(rankc = order(mcrps),
         rankw = order(mwis)) %>% 
  mutate(ranks = rankc == rankw, rankdist = abs(rankc - rankw)) %>% 
  arrange(rankw)
  

rank_cor <- matrix(.73, 
                   nrow = nrow(rank_scores), ncol = nrow(rank_scores))

for (i in 1:nrow(rank_scores)) {
  j <- rank_scores$rankc[i]
  rank_cor[i,j] <- rank_scores$corwc[j]
}

library(reshape2)

dat2 <- melt(rank_cor)
ggplot(dat2, aes(Var1, Var2, group=Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(fill = "CORR") +
  ylab("CRPS Rank") +
  xlab("WIS Rank") + 
  ylim(c(1,nrow(rank_scores))) +
  xlim(c(1,nrow(rank_scores))) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17))



ranked_scores %>% 
  ggplot() +
  geom_point(aes(x = rankw, y = rankc, colour = rankdist))


rank_scores %>% 
  filter(mcrps < 1) %>%
  ggplot() +
  geom_point(aes(x = mwis, y = mcrps, colour = rankdist), size = 2)





