setwd("~/quantile_fitting/")
source("./simulation/simulation_functions.R")
library(cmdstanr)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)

qgp_stan <- cmdstan_model(stan_file = 
                          './stan_models/cdf_quantile_normal_mix4.stan')

burn <- 6000
sample <- 7000

mod_loc <- "~/forecast-hub/FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- -1:3
# get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
# get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)




all_forecasts <- forecasts <- foreach(mod = models,
        .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                      "VGAM", "distr", "dplyr")
        ,.errorhandling = "remove"
        ,.combine = rbind) %:%
  foreach(date = sub_dates, .combine = rbind) %:%
    foreach(loc = locations) %:%
      foreach(h = horizons, .combine = rbind) %dopar% {
        source("./simulation_functions.R")
    
        if (dir.exists(paste0("model-fits/", mod)) == FALSE) {
          dir.create(paste0("model-fits/", mod))
        }
        
        if (dir.exists(paste0("model-fits/", mod, "/draws")) == FALSE) {
          dir.create(paste0("model-fits/", mod, "/draws"))
        }
        
        if (dir.exists(paste0("model-fits/", mod, "/forecasts")) == FALSE) {
          dir.create(paste0("model-fits/", mod, "/forecasts"))
        }
        
        if (dir.exists(paste0("model-fits/", mod, 
                              "/summary_diagnostics")) == FALSE) {
          dir.create(paste0("model-fits/", mod, "/summary_diagnostics"))
        }
        
        forc_file <- list.files(paste0(mod_loc, mod), pattern = date)
        all_forecasts <- read.csv(paste0(mod_loc, mod, "/", forc_file))
        
        forecast <- all_forecasts %>% 
          filter(location == as.character(loc), 
                 horizon == as.numeric(h), 
                 output_type == "quantile") %>% 
          unique()
        
        forecast %>%
          ggplot() +
          geom_point(aes(x = output_type_id, y = log(value + 1)))
        
        quantiles <- log(as.numeric(forecast$value) + 1)
        probs <- as.numeric(forecast$output_type_id)
        dat <- data.frame(quantile = quantiles,
                   prob = probs)
        
        dat <- dat %>% 
          filter(quantile != 0)
        
        stan_data <- make_stan_data(dat, size = 1, comps = 4)
        stan_samp <- qgp_stan$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            adapt_delta = .99,
                            refresh = 100)
        
        
        draws <- stan_samp$draws(format = "df")
        saveRDS(draws, paste0("model-fits/", mod, "/draws/", date, "-", mod,
                              ".rds"))
        
        summary <- stan_samp$summary()
        saveRDS(summary, paste0("model-fits/", mod, 
                                "/summary_diagnostics/", date, "-", mod, "_",
                              "summary.rds"))
        
        diagnostics <- stan_samp$diagnostic_summary()
        saveRDS(diagnostics, paste0("model-fits/", mod, 
                                "/summary_diagnostics/", date, "-", mod, "_",
                                "diagnostics.rds"))
        
        dist_samp <- draws$dist_samp
        
        est_quantiles <- quantile(dist_samp, probs = probs)
        est_quantiles <- ifelse(est_quantiles < 0, 0, est_quantiles)
        eval_quantiles <- ecdf(dist_samp)(quantiles)
        
        forecast$est_quantile <- exp(est_quantiles) - 1
        forecast$eval_quantile <- eval_quantiles
        
        forecast
  
    }

saveRDS(all_forecasts, "test_forecasts.rds")

