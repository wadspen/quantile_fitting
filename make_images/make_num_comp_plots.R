library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ka)

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
  summarise(muwd1 = mean(uwd1),
            muwd2 = mean(uwd2),
            mks = mean(ks),
            mkld = mean(kld),
            mtv = mean(tv)) %>% 
  mutate(comps = as.numeric(str_extract(model, "\\d+"))) %>% 
  mutate(n_known = str_detect(model, "n")) %>% 
  ungroup() %>% 
  select(-probs, -quants)


sum_scores %>% 
  filter(method == "MCMC", n == 1000, n_known == FALSE) %>% 
  ggplot() +
  geom_line(aes(x = comps, y = muwd1, colour = dist))# +
  # facet_wrap(~dist, scales = "free")



sum_scores %>% 
  filter(method == "variational", n == 1000, n_known == FALSE) %>% 
  ggplot() +
  geom_line(aes(x = comps, y = mtv, colour = dist))# +
# facet_wrap(~dist, scales = "free")


sum_scores %>% 
  filter(method == "MCMC", n == 1000, n_known == FALSE, dist == "gmix") %>% 
  arrange(comps) %>% 
  ungroup() %>% 
  mutate(eff = lag(muwd1)/muwd1)

sum_eff <- sum_scores %>% 
  arrange(method, n_known, dist, n, comps) %>% 
  group_by(method, n_known, dist, n) %>% 
  mutate(eff_uwd1 = lag(muwd1)/muwd1,
         eff_kld = lag(mkld)/mkld,
         eff_tv = lag(mtv)/mtv)


sum_eff %>% 
  filter(!is.na(eff_uwd1), method == "MCMC")





















