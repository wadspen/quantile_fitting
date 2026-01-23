setwd("../")
source("./simulation/simulation_functions.R")
library(dplyr)
library(scoringRules)
library(stringr)
library(evalcast)
library(lubridate)
library(evmix)
library(tidyr)
library(orderstats)
# library(StereoMorph)
library(cmdstanr)
library(dplyr)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
n.cores <- 64
n.cores <- 1
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)
print(paste("There are", n.cores, "cores!"))



args <- commandArgs()
forc_mod <- args[6]
fit_mod <- args[7]
h <- as.numeric(args[8])
horizon <- h
num_comps <- as.numeric(args[9])

hosp_loc <- "../FluSight-forecast-hub/target-data/target-hospital-admissions.csv"

hosp_data <- read.csv(hosp_loc) %>% 
  mutate(location = as.character(location)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0(0, location),
                           location)) %>%
  #select(-X) %>% 
  filter(year(date) >= 2023, date <= "2024-07-31")

mod_loc <- "./stan_models/"

mod_name <- ifelse(fit_mod == "clt_shs",
                   "cdf_quantile_normal_mixK_shs.stan",
                   ifelse(fit_mod == "clt", "cdf_quantile_normal_mixK.stan",
                          ifelse(fit_mod == "ord", "order_normal_mixK_quantiles.stan",
                                 ifelse(fit_mod == "clt_hs", 
                                        "cdf_quantile_normal_mixK_hs.stan",
                                        ifelse(fit_mod == "ord_hs",
                                               "order_normal_mixK_quantiles_hs.stan",
                                               
                                               ifelse(fit_mod == "ord_shs",
                                                      "order_normal_mixK_quantiles_shs.stan",
                                                      ifelse(fit_mod == "ind_shs",
                                                             "cdf_ind_quantile_normal_mixK_shs.stan",
                                                             ifelse(fit_mod == "clt_sb",
                                                                    "cdf_quantile_normal_mixK_sb.stan",
                                                                    ifelse(fit_mod == "ord_sb",
                                                                           "order_normal_mixK_quantiles_sb.stan",
                                                                           ifelse(fit_mod == "ind_sb", "cdf_ind_quantile_normal_mixK_sb.stan",
                                                                                  ifelse(fit_mod == "ind_hs", "cdf_ind_quantile_normal_mixK_hs.stan",
                                                                                         "cdf_ind_quantile_normal_mixK.stan")))))))))))

mod <- cmdstan_model(stan_file = paste0(mod_loc, mod_name))

burn <- 10000; burn <- 150
sample <- 50000;sample <- 100

mod_loc <- "../FluSight-forecast-hub/model-output/"
#models <- list.files(mod_loc)
#models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
sub_dates <- sub_dates[sub_dates <= "2024-07-31"]
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)
locations <- locations[locations != "US"]
locations <- c("US", locations)

#locations <- locations[3:5]
#sub_dates <- sub_dates[4:6]
#horizons <- 1:2
#models <- models[6:7]

date <- sub_dates[3]
loc <- locations[4]
# sub_dates <- "2023-11-18"
# locations <- "US"
# horizons <- 0
#all_forecasts <- foreach(date = sub_dates,
#                                      .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
#                                                    "VGAM", "distr", "dplyr",
#                                                    "lubridate", "scoringRules", "evalcast")
                                      # ,.errorhandling = "remove"
#                                      ,.combine = rbind) %:%
  #foreach(date = sub_dates, .combine = rbind) %:%
#  foreach(loc = locations, .combine = rbind) %dopar% {
  # foreach(h = horizons, .combine = rbind) %dopar% {
    source("./simulation/simulation_functions.R")
    #mod <- models[3]
    #h <- 1
    #loc <- locations[22]
    #date <- sub_dates[4]
 
    if (dir.exists(paste0("model-fits/", forc_mod, "/results")) == FALSE) {
      dir.create(paste0("model-fits/", forc_mod, "/results"), recursive = TRUE)
    }
    
    
    if (dir.exists(paste0("model-fits/", forc_mod, 
                          "/summary_diagnostics")) == FALSE) {
      dir.create(paste0("model-fits/", forc_mod, "/summary_diagnostics"),
                 recursive = TRUE)
    }
   print(paste0(mod_loc, forc_mod)) 
    forc_file <- list.files(paste0(mod_loc, forc_mod), pattern = date)
    forecasts <- read.csv(paste0(mod_loc, forc_mod, "/", forc_file))
 print("does it get here?")      
    forecast <- forecasts %>% 
      filter(location == as.character(loc) | location == as.numeric(loc), 
             horizon == as.numeric(h) | horizon == as.character(h), 
             output_type == "quantile") %>% 
      unique()
    
    #print(head(forecast)) 
    quantiles <- log(as.numeric(forecast$value) + 1)
    probs <- as.numeric(forecast$output_type_id)
    dat <- data.frame(quantile = quantiles,
                      prob = probs)
    
    dat <- dat %>% 
      filter(quantile != 0)
    print("is it broken here")
    stan_data <- make_stan_data(dat, size = 1, comps = num_comps, alph = 1)
    
    if (loc == "US" & date == "2023-11-18" & as.numeric(h) == 0) {
      
      stan_samp <- mod$sample(data = stan_data,
                              iter_warmup = burn,
                              iter_sampling = sample,
                              chains = 4,
                              adapt_delta = .999,
                              refresh = 100)    
      draws <- stan_samp$draws(format = "df") %>% 
        filter(if_all(everything(), ~ !is.na(.)))
      
      summary <- stan_samp$summary() %>% filter(variable == "dist_samp")
      diag_list <- lapply(stan_samp$diagnostic_summary(), FUN = mean)
      diag <- diag_list %>% data.frame()
      diag$max_rhat <- max(summary$rhat, na.rm = TRUE)
      diag$min_ess <- min(summary$ess_bulk, na.rm = TRUE)
      diag$min_esst <- min(summary$ess_tail, na.rm = TRUE)
      
      
      saveRDS(diag, paste0("model-fits/", forc_mod, 
                           "/summary_diagnostics/diagnostics.rds"))
      
      # saveRDS(diag, paste0("model-fits-ord/all_diags/", mod, "diagnostics.rds"))
      
      draws <- stan_samp$draws(format = "df") %>% 
        filter(if_all(everything(), ~ !is.na(.)))
      
      drawsQ <- draws %>% 
        dplyr::select(contains("Q_rep"))
      
      
      data_sum <- dat %>%
        mutate(date = date) %>% 
        left_join(hosp_data %>% 
                    filter(location == loc),
                  by = "date") %>% 
        mutate(reference_date = date) 
      
      data_sum <- data_sum %>% 
        dplyr::select(-date) %>% 
        mutate(true_value = log(value + 1),
               date = date(reference_date) + h * 7) %>% 
        mutate(low95 = apply(drawsQ, MARGIN = 2, 
                             FUN = quantile, probs = 0.025),
               upp95 = apply(drawsQ, MARGIN = 2, 
                             FUN = quantile, probs = 0.975),
               mean = apply(drawsQ, MARGIN = 2,
                            FUN = mean),
               mae = abs(quantile - mean),
               mse = (quantile - mean)^2
        ) %>% 
        mutate(cover = between(quantile, low95, upp95))
      
      data_sum$est_quantile <- quantile(draws$dist_samp, probs = data_sum$prob)
      
      score_sum <- data_sum %>% 
        group_by(location, reference_date, date) %>% 
        summarise(wis = weighted_interval_score(prob, quantile, 
                                                unique(true_value)),
                  est_wis = weighted_interval_score(prob, est_quantile,
                                                    unique(true_value)),
                  mae = mean(mae),
                  mse = mean(mse),
                  coverm = mean(cover))
      
      
      score_sum$crps <- crps_sample(unique(data_sum$true_value), draws$dist_samp)
      score_sum$logs <- logs_sample(unique(data_sum$true_value), draws$dist_samp)
      score_sum <- score_sum %>% 
        mutate(horizon = h,
               forc_model = forc_mod,
               fit_model = fit_mod)
      
      score_sum
    } else {
      stan_samp <- mod$sample(data = stan_data,
                              iter_warmup = burn,
                              iter_sampling = sample,
                              chains = 1,
                              adapt_delta = .999,
                              refresh = 100)
      
      
      draws <- stan_samp$draws(format = "df") %>% 
        filter(if_all(everything(), ~ !is.na(.)))
      
      drawsQ <- draws %>% 
        dplyr::select(contains("Q_rep"))
      
      
      data_sum <- dat %>%
        mutate(date = date) %>% 
        left_join(hosp_data %>% 
                    filter(location == loc),
                  by = "date") %>% 
        mutate(reference_date = date) 
     
      data_sum <- data_sum %>% 
        dplyr::select(-date) %>% 
        mutate(true_value = log(value + 1),
               date = date(reference_date) + h * 7) %>% 
        mutate(low95 = apply(drawsQ, MARGIN = 2, 
                             FUN = quantile, probs = 0.025),
               upp95 = apply(drawsQ, MARGIN = 2, 
                             FUN = quantile, probs = 0.975),
               mean = apply(drawsQ, MARGIN = 2,
                            FUN = mean),
               mae = abs(quantile - mean),
               mse = (quantile - mean)^2
        ) %>% 
        mutate(cover = between(quantile, low95, upp95))
      
      data_sum$est_quantile <- quantile(draws$dist_samp, probs = data_sum$prob)
      
      score_sum <- data_sum %>% 
        group_by(location, reference_date, date) %>% 
        summarise(wis = weighted_interval_score(prob, quantile, 
                                                unique(true_value)),
                  est_wis = weighted_interval_score(prob, est_quantile,
                                                    unique(true_value)),
                  mae = mean(mae),
                  mse = mean(mse),
                  coverm = mean(cover))
      
      
      score_sum$crps <- crps_sample(unique(data_sum$true_value), draws$dist_samp)
      score_sum$logs <- logs_sample(unique(data_sum$true_value), draws$dist_samp)
      score_sum <- score_sum %>% 
        mutate(horizon = h,
               forc_model = forc_mod,
               fit_model = fit_mod)
      
      score_sum
      
      
    }
    
    
  #}

saveRDS(all_forecasts, paste0("./model-fits/", forc_mod, "/results/horizon_",
                              h, "_", fit_mod, "_ncomp_", num_comps, ".rds"))





