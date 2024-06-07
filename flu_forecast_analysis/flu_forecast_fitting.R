source("~/flu_research/nonlinear_flu_forecast/get_data_functions.R")
setwd("~/quantile_fitting/")
library(cmdstanr)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)

qgp_stan <- cmdstan_model(stan_file = 
                          './stan_models/normal_t_mix4_quantiles.stan')

burn <- 6000
sample <- 7000

mod_loc <- "~/forecast-hub/FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)




forecasts <- foreach(mod = models,
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
        # dat <- dat[-c(21:23),]
        
        stan_data <- make_stan_data(dat, size = 1, comps = 5)
        stan_samps <- qgp_stan$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            adapt_delta = .9999,
                            refresh = 100)
        
        
        draws <- stan_samps$draws(format = "df")
        saveRDS(draws, paste0("model-fits/", mod, "/draws/", date, "-", mod,
                              ".rds"))
        dist_samp <- draws$dist_samps
        
        
        # dist_samp <- ifelse(dist_samp < 0, 0, dist_samp)
        est_quantiles <- quantile(dist_samp, probs = probs)
        est_quantiles <- ifelse(est_quantiles < 0, 0, est_quantiles)
        # q <- q[-length(q)]
        # plot(q ~ probs)#[-length(probs)])
        # plot(q[-c(21:23)] ~ stan_data$Q)
        # abline(a = 0, b = 1)
        eval_quantiles <- ecdf(dist_samp)(quantiles)
        
        forecast$est_quantile <- exp(est_quantiles) - 1
        forecast$eval_quantile <- eval_quantiles
        
        forecast
        # plot(test~probs)
        # abline(a = 0, b = 1)
        # sum((test[6:22] - probs[6:22])^2) 
                              #Flu-base, US, date 2
                              #.033299 no 0
                              #.03203 keep only smallest 0
                              #.070431 keep all 0s
                              #.085867 keep only largest 0
        
                              #"CEPH-Rtrend_fluH" 17 date2
                              #all data .033506
                              #remove all 0s .009452
                              #keep only smallest 0 .01302
                              #keep largest 0
        
                              #GT-FluFNP 17 2
                              #keep all data .252011
                              #remove all 0s .171061
                              #keep only largest 0 .127768
                              #keep only smallest 0 .071998
    
    
      }



