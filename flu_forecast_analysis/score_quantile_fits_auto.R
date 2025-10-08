source("../simulation/simulation_functions.R")
source("./get_fit_cover.R")
library(dplyr)
library(scoringRules)
library(stringr)
library(evalcast)
library(lubridate)
library(evmix)
library(tidyr)
library(orderstats)
library(StereoMorph)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
n.cores <- 64
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)


args <- commandArgs()
mod <- args[6]
method <- args[7]
print(mod)
print(method)
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
# mod_loc <- "../model-fits-ord/"
mod_loc <- ifelse(str_detect(method, "ind"), "../model-fits-ind/",
                  ifelse(str_detect(method, "ord"), "../model-fits-ord/",
                                    "../model-fits/"))

indt <- ifelse(str_detect(mod_loc, "ind"), TRUE, FALSE)
ordt <- ifelse(str_detect(mod_loc, "ord"), TRUE, FALSE)

models <- list.files(forc_loc)
sub_dates <- substr(list.files(paste0(forc_loc, 
                                "FluSight-baseline/")),
                    1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(forc_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(forc_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)
print(locations)
hosp_data <- read.csv(hosp_loc) %>% 
  mutate(location = as.character(location)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0(0, location),
			   location)) %>%
  #select(-X) %>% 
  filter(year(date) >= 2023)

# mod <- "CEPH-Rtrend_fluH"
# h <- 2
# loc <- "US"
# sub_date <- "2023-10-14"

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
	  mutate(location = as.character(location)) %>%
  	  mutate(location = ifelse(nchar(location) < 2, paste0(0, location),
			   location)) %>%	     
          filter(location == loc, 
                 horizon == as.numeric(h), 
                 output_type == "quantile") %>% 
          unique()

  	forecast <- forecast %>%
		mutate(date = date(reference_date) + 7*as.numeric(horizon),
		       location = location)
        
        quantiles <- log(as.numeric(forecast$value) + 1)
        probs <- as.numeric(forecast$output_type_id)
    
        draws <- readRDS(paste0(mod_loc, mod, "/draws/", sub_date, "-", 
                                loc, "-", h, "-", mod, ".rds"))
        
        # qcorr <- make_qcorr(probs)
        # tot_coverage <- get_fit_coverage(draws, 
        #                  true_quantiles = quantiles, true_probs = probs,
        #                  n = NULL, n_known = FALSE, n_modeled = TRUE,
        #                  ind = indt, order = ordt, QCorr = qcorr, 
        #                  num_samps = 400)
        # 
        # cover <- apply(tot_coverage, MARGIN = 2, FUN = mean) %>% 
        #   t() %>% 
        #   as.data.frame()
        #cover <- readRDS(paste0(mod_loc, mod, "/coverage/", sub_date, "-", 
        #                        loc, "-", h, "-", mod, ".rds"))
        #cover <- cover %>%
	#	dplyr::select(contains("cover"), contains("wid"), prob) %>%
	#	filter(as.numeric(prob) %in% c(0.01, 0.99, 0.025, 0.975, 0.25, 0.75)) %>%
	#	dplyr::select(contains("50"), contains("95"), contains("98"), prob) %>%
	#	pivot_wider(values_from = 1:6, names_from = prob)


  true_hosp <- hosp_data %>% 
    # dplyr::select(-X) %>%
    filter(date == date(sub_date) + 7*h, location == loc) %>% 
    mutate(value = log(value + 1), date = date(date)) %>% 
    mutate(true_value = value) %>%
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
			ssq = sum((log(value + 1) - log(est_quantile + 1))^2),
			sap = sum(abs(eval_quantile - output_type_id)),
			saq = sum(abs(log(value + 1) - log(est_quantile + 1))),
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
			ssq = sum((log(value + 1) - log(est_quantile + 1))^2),
			sap = sum(abs(eval_quantile - output_type_id)),
			saq = sum(abs(log(value + 1) - log(est_quantile + 1))),
			int = get_lm(est_quantile, value, stat = "int"),
			se_int = get_lm(est_quantile, value, stat = "se_int"),
			slope = get_lm(est_quantile, value, stat = "se_slope"),
			se_slope = get_lm(est_quantile, value, stat = "df"),
			n = dplyr::n()
		)
	
	all_score <- data.frame(model = mod, date = sub_date, location = loc, horizon = h, 
			     logs = logs, crps = crps, scores, scores0)
        
	all_score    
         
}

        saveRDS(all_scores, paste0(mod_loc, mod, 
                                "/scores/all_scores_no_cover.rds"))
        























