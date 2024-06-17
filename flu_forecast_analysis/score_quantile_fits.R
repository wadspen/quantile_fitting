source("../simulation/simulation_functions.R")
library(dplyr)
library(scoringRules)
library(stringr)
library(evalcast)
library(lubridate)
library(evmix)
library(tidyr)
library(StereoMorph)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
#my.cluster <- makeCluster(n.cores, type = "PSOCK")
#doParallel::registerDoParallel(cl = my.cluster)
#foreach::getDoParRegistered()
#foreach::getDoParWorkers()
#registerDoMC(cores = n.cores)
#

args <- commandArgs()
mod <- args[6]
print(mod)
eval_dist <- function(ests, probs, p = 1) {
	ests <- c(ests, 1)
	probs <- c(probs, 1)
	for_diffs <- c(0, probs)

	#est_diffs <- distancePointToPoint(ests)
	prob_diffs <- distancePointToPoint(for_diffs)
	return((p + 1)*sum(abs(ests*prob_diffs - probs*prob_diffs)^p)^(1/p))
}

get_lm <- function(y, x, stat = "slope") {
	mod <- lm(y ~ x)
	if (stat == "int") {
		return(mod$coefficient[1])
	} else if (stat == "slope") {
		return(mod$coefficient[2])
	} else if (stat == "df") {
		return(mod$df)
	} else if (stat == "se_int") {
		return(coef(summary(mod))[, "Std. Error"][1])
	} else if (stat == "se_slope") {
		return(coef(summary(mod))[, "Std. Error"][2])
	}
}



forc_loc <- "../../../FluSight-forecast-hub/model-output/"
hosp_loc <- paste0(forc_loc, "../target-data/target-hospital-admissions.csv")
mod_loc <- "../model-fits/"
models <- list.files(forc_loc)
sub_dates <- substr(list.files(paste0(forc_loc, 
                                "FluSight-baseline/")),
                    1, 10)
horizons <- -1:3
get_loc_file <- list.files(paste0(forc_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(forc_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

hosp_data <- read.csv(hosp_loc) %>% 
  #select(-X) %>% 
  filter(year(date) >= 2023)
#
#mod <- "CEPH-Rtrend_fluH"
#h <- 2
#loc <- "17"
#sub_date <- "2023-10-14"

all_scores <- foreach(sub_date = sub_dates,
                                       .packages = c("distr", "dplyr", 
                                                     "stringr", "scoringRules",
                                                     "evalcast", "StereoMorph", "tidyr")
                                       ,.errorhandling = "remove"
                                       ,.combine = rbind) %:%
   #foreach(sub_date = sub_dates, .combine = rbind) %:%
     foreach(loc = locations, .combine = rbind) %:%
       foreach(h = horizons, .combine = rbind) %dopar% {
        
       if (dir.exists(paste0(mod_loc, mod, "/scores")) == FALSE) {
          dir.create(paste0(mod_loc, mod, "/scores"))
        }

        forc_file <- paste0(mod_loc, mod, "/forecasts/", sub_date, "-",
			    loc, "-", h, "-", mod, ".rds")

	forcRDS <- readRDS(forc_file)
        forecast <- forcRDS %>% 
          filter(location == as.numeric(loc), 
                 horizon == as.numeric(h), 
                 output_type == "quantile") %>% 
          unique()

  	forecast <- forecast %>%
		mutate(date = date(reference_date) + 7*as.numeric(horizon),
		       location = as.numeric(location))
        
        quantiles <- log(as.numeric(forecast$value) + 1)
        probs <- as.numeric(forecast$output_type_id)
    
        draws <- readRDS(paste0(mod_loc, mod, "/draws/", sub_date, "-", 
                                loc, "-", h, "-", mod, ".rds"))
       
        cover <- readRDS(paste0(mod_loc, mod, "/coverage/", sub_date, "-", 
                                loc, "-", h, "-", mod, ".rds"))
        cover <- cover %>%
		dplyr::select(contains("cover"), contains("wid"), prob) %>%
		filter(as.numeric(prob) %in% c(0.01, 0.99, 0.025, 0.975, 0.25, 0.75)) %>%
		dplyr::select(contains("50"), contains("95"), contains("98"), prob) %>%
		pivot_wider(values_from = 1:6, names_from = prob)


        true_hosp <- hosp_data %>% 
	  dplyr::select(-X) %>%
          filter(date == date(sub_date) + 7*h, location == loc) %>% 
          mutate(value = log(value + 1), date = date(date),
	         location = as.numeric(location), true_value = value) %>%
	  dplyr::select(-value)
        
        crps <- crps_sample(true_hosp$true_value, draws$dist_samp)
        logs <- logs_sample(true_hosp$true_value, draws$dist_samp)
        
	forecast <- forecast %>%
		left_join(true_hosp, by = c("date", "location")) %>%
		mutate(output_type_id = as.numeric(output_type_id),
		       value = as.numeric(value))


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
			n = n()
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
			n = n()
		)
	
	scores <- data.frame(model = mod, date = sub_date, location = loc, horizon = h, 
			     logs = logs, crps = crps, scores, scores0, cover)
        
	scores    
     
       }

        saveRDS(all_scores, paste0(mod_loc, mod, 
                                "/scores/all_scores.rds"))
        























