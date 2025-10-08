setwd("../")
source("./simulation/simulation_functions.R")
library(cmdstanr)
library(dplyr)
library(evmix)
library(distfromq)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
n.cores <- 64
#n.cores <- 1
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)
print(paste("There are", n.cores, "cores!"))

args <- commandArgs()
mod <- args[6]


mod_loc <- "../FluSight-forecast-hub/model-output/"
#models <- list.files(mod_loc)
#models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

forc_loc <- "../../FluSight-forecast-hub/model-output/"
hosp_loc <- paste0(forc_loc, "../target-data/target-hospital-admissions.csv")

hosp_loc <- paste0(forc_loc, "../target-data/target-hospital-admissions.csv")
hosp_data <- read.csv(hosp_loc) %>% 
  mutate(location = as.character(location)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0(0, location),
                           location)) %>%
  #select(-X) %>% 
  filter(year(date) >= 2023)

#locations <- locations[3:5]
#sub_dates <- sub_dates[4:6]
#horizons <- 1:2
#models <- models[6:7]


#sub_dates <- "2023-11-18"
#locations <- "US"
#horizons <- 0
n <- 60000
all_forecasts <- forecasts <- foreach(date = sub_dates,
        .packages = c("evmix", "distfromq", "EnvStats",
                      "VGAM", "distr", "dplyr")
        ,.errorhandling = "remove"
        ,.combine = rbind) %:%
  #foreach(date = sub_dates, .combine = rbind) %:%
    foreach(loc = locations, .combine = rbind) %:%
      foreach(h = horizons, .combine = rbind) %dopar% {
        source("./simulation/simulation_functions.R")
   #mod <- models[3]
   #h <- 1
   #loc <- locations[22]
   #date <- sub_dates[4]
        if (dir.exists(paste0("model-fits-ord/", mod)) == FALSE) {
          dir.create(paste0("model-fits-ord/", mod))
        }
        
        if (dir.exists(paste0("model-fits-ord/", mod, "/draws")) == FALSE) {
          dir.create(paste0("model-fits-ord/", mod, "/draws"))
        }
        
        if (dir.exists(paste0("model-fits-ord/", mod, "/forecasts")) == FALSE) {
          dir.create(paste0("model-fits-ord/", mod, "/forecasts"))
        }
        
        if (dir.exists(paste0("model-fits-ord/", mod, 
                              "/summary_diagnostics")) == FALSE) {
          dir.create(paste0("model-fits-ord/", mod, "/summary_diagnostics"))
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
        
        
        true_hosp <- hosp_data %>% 
          # dplyr::select(-X) %>%
          filter(date == date(sub_date) + 7*h, location == loc) %>% 
          mutate(value = log(value + 1), date = date(date),
                 true_value = value) %>%
          dplyr::select(-value)
    
        qspline <- make_q_fn(probs, quantiles)
        rspline <- make_r_fn(probs, quantiles)
        
        sp_samp <- rspline(n)
        
        crps_spline <- crps_sample(true_hosp$true_value, sp_samp)
        logs_spline <- logs_sample(true_hosp$true_value, sp_samp)
        
        qkern <- function(p) {qkden(p, quantiles, kernel = "gaussian")}
        rkern <- function(n) {rkden(n, quantiles, kernel = "gaussian")}
        
        k_samp <- rkern(n)
        
        crps_kern <- crps_sample(true_hosp$true_value, k_samp)
        logs_kern <- logs_sample(true_hosp$true_value, k_samp)
        
        
        
        
        scores0 <- forecast %>%
          filter(value != 0) %>%
          ungroup() %>%
          mutate(output_type_id = as.numeric(output_type_id)) %>%
          reframe(
            true_wis = weighted_interval_score(as.numeric(output_type_id),
                                               log(as.numeric(value) + 1),
                                               as.numeric(unique(true_value))),
            est_wis = weighted_interval_score(as.numeric(output_type_id),
                                              log(as.numeric(est_quantile) + 1),
                                              as.numeric(unique(true_value))),
            dist_diff = eval_dist(eval_quantile, as.numeric(output_type_id)),
            dist_diff2 = eval_dist(eval_quantile, as.numeric(output_type_id), p = 2),
            ssp = sum((eval_quantile - output_type_id)^2),
            ssq = sum((value - est_quantile)^2),
            sap = sum(abs(eval_quantile - output_type_id)),
            saq = sum(abs(value - est_quantile)),
            int = get_lm(est_quantile, value, stat = "int"),
            se_int = get_lm(est_quantile, value, stat = "se_int"),
            slope = get_lm(est_quantile, value, stat = "se_slope"),
            se_slope = get_lm(est_quantile, value, stat = "df"),
            n = dplyr::n()
          )	
        colnames(scores0) <- paste0(colnames(scores0), 0)
        
        scores <- forecast %>%
          ungroup() %>%
          reframe(
            true_wis = weighted_interval_score(as.numeric(output_type_id),
                                               log(as.numeric(value) + 1),
                                               as.numeric(unique(true_value))),
            est_wis = weighted_interval_score(as.numeric(output_type_id),
                                              log(as.numeric(est_quantile) + 1),
                                              as.numeric(unique(true_value))),
            dist_diff = eval_dist(eval_quantile, output_type_id),
            dist_diff2 = eval_dist(eval_quantile, output_type_id, p = 2),
            ssp = sum((eval_quantile - output_type_id)^2),
            ssq = sum((value - est_quantile)^2),
            sap = sum(abs(eval_quantile - output_type_id)),
            saq = sum(abs(value - est_quantile)),
            int = get_lm(est_quantile, value, stat = "int"),
            se_int = get_lm(est_quantile, value, stat = "se_int"),
            slope = get_lm(est_quantile, value, stat = "se_slope"),
            se_slope = get_lm(est_quantile, value, stat = "df"),
            n = dplyr::n()
          )
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
	saveRDS(draws, paste0("model-fits-duh/", mod, "/draws/", date, "-", loc,
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
        saveRDS(forecast, paste0("model-fits-ord/", mod,
				 "/forecasts/", date, "-", loc,
				 "-", h, "-", mod, ".rds"))

  
    }




