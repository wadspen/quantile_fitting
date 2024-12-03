setwd("../")
source("./simulation/simulation_functions.R")
library(cmdstanr)
library(dplyr)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
n.cores <- 1
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)
print(paste("There are", n.cores, "cores!"))

args <- commandArgs()
mod <- args[6]
print(mod)
qgp_stan <- cmdstan_model(stan_file = 
                          './stan_models/cdf_quantile_normal_mix4.stan')

burn <- 20000; # burn <- 150
sample <- 60000;# sample <- 100

mod_loc <- "../FluSight-forecast-hub/model-output/"
#models <- list.files(mod_loc)
#models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

#locations <- locations[3:5]
#sub_dates <- sub_dates[4:6]
#horizons <- 1:2
#models <- models[6:7]

locations = "US"
sub_dates = "2023-11-18"
burn <- 1500
sample = 1500
loc = "US"
date = "2023-11-18"
h = 0
#all_forecasts <- forecasts <- foreach(date = sub_dates,
#        .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
#                      "VGAM", "distr", "dplyr")
#        #,.errorhandling = "remove"
#        ,.combine = rbind) %:%
#  #foreach(date = sub_dates, .combine = rbind) %:%
#    foreach(loc = locations, .combine = rbind) %:%
#      foreach(h = horizons, .combine = rbind) %dopar% {
        source("./simulation/simulation_functions.R")
   #mod <- models[3]
   #h <- 1
   #loc <- locations[22]
   #date <- sub_dates[4]
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
        forecasts <- read.csv(paste0(mod_loc, mod, "/", forc_file))
        
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
    
        stan_data <- make_stan_data(dat, size = 1, comps = 4)
        if (loc == "US" & date == "2023-11-18" & as.numeric(h) == 0) {
	   
	    stan_samp <- qgp_stan$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 4,
                            adapt_delta = .999,
                            refresh = 100)    
	    draws <- stan_samp$draws(format = "df")

	    summary <- stan_samp$summary() %>% filter(variable == "dist_samp")
	    diag_list <- lapply(stan_samp$diagnostic_summary(), FUN = mean)
	    diag <- diag_list %>% data.frame()
	    diag$max_rhat <- max(summary$rhat, na.rm = TRUE)
	    diag$min_ess <- min(summary$ess_bulk, na.rm = TRUE)
	    diag$min_esst <- min(summary$ess_tail, na.rm = TRUE)
	    
	    saveRDS(diag, paste0("model-fits/all_diags/", mod, "diagnostics.rds"))
    } else {
        stan_samp <- qgp_stan$sample(data = stan_data,
                            iter_warmup = burn,
                            iter_sampling = sample,
                            chains = 1,
                            adapt_delta = .999,
                            refresh = 100)
        
        
        draws <- stan_samp$draws(format = "df")
     }
	saveRDS(draws, paste0("model-fits/", mod, "/draws/", date, "-", loc,
			      "-", h, "-", mod,
                              ".rds"))
        
        #summary <- stan_samp$summary()
        #saveRDS(summary, paste0("model-fits/", mod, 
        #                        "/summary_diagnostics/", date, "-", loc,
	#			"-", h, "-", mod, "_",
        #                      "summary.rds"))
        #
        #diagnostics <- stan_samp$diagnostic_summary()
        #saveRDS(diagnostics, paste0("model-fits/", mod, 
        #                        "/summary_diagnostics/", date, "-", loc,
	#			"-", h, "-", mod, "_",
        #                        "diagnostics.rds"))
        
        dist_samp <- draws$dist_samp
        
        est_quantiles <- quantile(dist_samp, probs = probs)
        est_quantiles <- ifelse(est_quantiles < 0, 0, est_quantiles)
        eval_quantiles <- ecdf(dist_samp)(quantiles)
        
        forecast$est_quantile <- exp(est_quantiles) - 1
        forecast$eval_quantile <- eval_quantiles
        forecast$model <- mod 
        saveRDS(forecast, paste0("model-fits/", mod,
				 "/forecasts/", date, "-", loc,
				 "-", h, "-", mod, ".rds"))

  
#    }
#write.csv(all_forecasts, "test_forecasts.csv") 
#saveRDS(all_forecasts, "test_forecasts.rds")

#forecasts <- all_forecasts %>%
#	group_split(model, reference_date)
#
#
#for (i in 1:length(forecasts)) {
#        forc <- forecasts[[i]]
#	ref_date <- unique(forc$reference_date)
#	model <- unique(forc$model)
#	forc <- forc %>%
#		dplyr::select(-model)
#	
#	save_name <- paste0("./model-fits/", model, "/forecasts/",
#	       ref_date, "-", model, ".csv")
#
#	write.csv(forc, save_name, row.names = FALSE)
#}






