setwd("~/quantile_fitting/")
source("./simulation/simulation_functions.R")
library(dplyr)
library(lubridate)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)



read.csv("~/forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv")
forc_loc <- "../../FluSight-forecast-hub/model-output/"
forc_loc <- "~/forecast-hub/FluSight-forecast-hub/model-output/"
hosp_loc <- paste0(forc_loc, "../target-data/target-hospital-admissions.csv")
mod_loc <- "./model-fits/"
mod_loc <- "~/forecast-hub/FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
sub_dates <- substr(list.files(paste0(mod_loc, 
                                "FluSight-forecast-hub/model-output/")),
                    1, 10)
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- -1:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

hosp_data <- read.csv(hosp_loc) %>% 
  select(-X) %>% 
  filter(year(date) >= 2023)

mod <- "FluSight-lop_norm"
h <- 2
loc <- "17"
sub_date <- "2023-10-28"

# all_forecasts <- forecasts <- foreach(mod = models,
#                                       .packages = c("distr", "dplyr", 
#                                                     "stringr", "scoringRules",
#                                                     "eval_cast")
#                                       ,.errorhandling = "remove"
#                                       ,.combine = rbind) %:%
#   foreach(sub_date = sub_dates, .combine = rbind) %:%
#     foreach(loc = locations) %:%
#       foreach(h = horizons, .combine = rbind) %dopar% {
        
        forc_file <- list.files(paste0(forc_loc, mod), pattern = sub_date)
        all_forecasts <- read.csv(paste0(forc_loc, mod, "/", forc_file))
        
        forecast <- all_forecasts %>% 
          filter(location == as.character(loc), 
                 horizon == as.numeric(h), 
                 output_type == "quantile") %>% 
          unique()
        
        quantiles <- log(as.numeric(forecast$value) + 1)
        probs <- as.numeric(forecast$output_type_id)
    
        draws <- readRDS(paste0(mod_loc, mod, "/draws/", sub_date, "-", 
                                loc, "-", h, "-", mod, ".rds"))
        
        true_hosp <- hosp_data %>% 
          filter(date == date(sub_date) + 7*h, location == loc) %>% 
          mutate(value = log(value + 1))
        
        crps <- crps_sample(true_hosp$value, draws$dist_samp)
        logs <- logs_sample(true_hosp$value, draws$dist_samp)
        
        
        all_pis <- draws[, str_detect(colnames(draws), "pi")]
        all_mus <- draws[, str_detect(colnames(draws), "mus")]
        all_sigmas <- draws[, str_detect(colnames(draws), "sigmas")]
        all_ns <- draws$n
        
        num_samps <- nrow(draws)
        
        m <- 1
        samp_quantiles <- matrix(NA, nrow = num_samps, ncol = length(probs))
        repeat{
          num <- sample(num_samps, 1)
          mus <- unlist(all_mus[num,])
          sigmas <- unlist(all_sigmas[num,])
          pi <- unlist(all_pis[num,])
          pi[which(pi < 0)] <- 0
          if (sum(pi) < 1) {
            pi[which.max(pi)] <- pi[which.max(pi)] + (1 - sum(pi))
          } else if (sum(pi) > 1) {
            pi[which.min(pi)] <- pi[which.min(pi)] + (1 - sum(pi))
          }
          n <- unlist(all_ns[num])
          
          normmix <- UnivarMixingDistribution(Norm(mus[1], sigmas[1]),
                                              Norm(mus[2], sigmas[2]), 
                                              Norm(mus[3], sigmas[3]), 
                                              Norm(mus[4], sigmas[4]), 
                                              mixCoeff = pi)
          samp <- r(normmix)(n)
          samp_quantiles[m,] <- quantile(samp, probs = probs)
          
          m <- m + 1
          if (m > num_samps) {break}
        }
        
        
        quant_bounds <- apply(samp_quantiles, MARGIN = 2, 
                              FUN = quantile, 
                              probs = probs, na.rm = TRUE)
        
        quant_bounds <- data.frame(t(quant_bounds))
        colnames(quant_bounds) <- as.character(probs)
        quant_bounds$prob <- probs
        quant_bounds$quantile <- quantiles
        
        quant_bounds <- quant_bounds %>% 
          mutate(cover98 = between(quantile, `0.01`, `0.99`),
                 cover95 = between(quantile, `0.025`, `0.975`), 
                 cover90 = between(quantile, `0.05`, `0.95`),
                 cover80 = between(quantile, `0.1`, `0.9`),
                 cover70 = between(quantile, `0.15`, `0.85`),
                 cover60 = between(quantile, `0.2`, `0.8`),
                 cover50 = between(quantile, `0.25`, `0.75`),
                 cover40 = between(quantile, `0.3`, `0.7`),
                 cover30 = between(quantile, `0.35`, `0.65`),
                 cover20 = between(quantile, `0.4`, `0.6`),
                 cover10 = between(quantile, `0.45`, `0.55`))
    
        
        # quant_bounds %>% 
        #   filter(between(prob, .25, .75)) %>%
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
        #   theme_bw()
    
    
    
      # }
























