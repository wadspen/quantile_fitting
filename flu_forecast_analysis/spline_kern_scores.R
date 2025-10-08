setwd("../")
source("./simulation/simulation_functions.R")
source("./flu_forecast_analysis/get_fit_cover.R")
library(cmdstanr)
library(dplyr)
library(evmix)
library(lubridate)
library(distfromq)
library(VGAM)
library(EnvStats)
library(distr)
library(evalcast)
library(scoringRules)
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
print(paste("There are", n.cores, "cores!"))

args <- commandArgs()
mod <- args[6]
print(mod)

forc_loc <- "./../../FluSight-forecast-hub/model-output/"
hosp_loc <- paste0(forc_loc, "../target-data/target-hospital-admissions.csv")

hosp_data <- read.csv(hosp_loc) %>% 
  mutate(location = as.character(location)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0(0, location),
                           location)) %>%
  #select(-X) %>% 
  filter(year(date) >= 2023)


mod_loc <- "../FluSight-forecast-hub/model-output/"
models <- list.files(mod_loc)
models <- models[models != "README.md"]
sub_dates <- substr(list.files(paste0(mod_loc, "FluSight-baseline")), 1, 10)
print(sub_dates)
horizons <- 0:3
get_loc_file <- list.files(paste0(mod_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(mod_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

#locations <- locations[3:5]
#sub_dates <- sub_dates[4:6]
#horizons <- 1:2
#models <- models[6:7]


#sub_dates <- "2023-11-18"
#locations <- "US"
#horizons <- 0
lsd <- length(sub_dates)
print(sub_dates)
# sub_dates <- sub_dates[(lsd - 5):lsd]
print(sub_dates)
sample_size <- 5000
all_scores <- forecasts <- foreach(sub_date = sub_dates,
                                      .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                                    "VGAM", "distr", "dplyr",
                                                    "scoringRules", "evalcast",
                                                    "StereoMorph")
                                      ,.errorhandling = "remove"
                                      ,.combine = rbind) %:%
  foreach(loc = locations, .combine = rbind) %:%
  # foreach(mod = models, .combine = rbind) %:%
  foreach(h = horizons, .combine = rbind) %dopar% {
    source("./simulation/simulation_functions.R")
    # mod <- models[3]
    # h <- 1
    # loc <- locations[22]
    # sub_date <- sub_dates[4]
    true_hosp <- hosp_data %>% 
      # dplyr::select(-X) %>%
      filter(date == date(sub_date) + 7*h, location == loc) %>% 
      mutate(value = log(value + 1), date = date(date)) %>% 
      mutate(true_value = value) %>%
      dplyr::select(-value)
    
    forc_file <- list.files(paste0(mod_loc, mod), pattern = sub_date)
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
    
    
    forecast0 <- forecast %>% 
      filter(value != 0)
    
    dat <- dat %>% 
      filter(quantile != 0)
    
    probs0 <- dat$prob
    quantiles0 <- dat$quantile
    
    qspline <- make_q_fn(probs0, quantiles0)
    pspline <- make_p_fn(probs0, quantiles0)
    puspline <- function(x) {pdist(qspline(x))}
    rspline <- make_r_fn(probs0, quantiles0)
    
    eval_quantile_spline <- pspline(quantiles0)
    
    spline_draws <- rspline(sample_size)
    spline_quant <- qspline(probs0)
    print("does it get here?")
    spline_scores0 <- data.frame(
      model = mod,
      date = sub_date, 
      location = loc, 
      horizon = h, 
      method = "spline",
      zeros = TRUE,
      crps = scoringRules::crps_sample(true_hosp$true_value, spline_draws),
      logs = scoringRules::logs_sample(true_hosp$true_value, spline_draws),
      
      true_wis = weighted_interval_score(as.numeric(probs0),
                                         quantiles0,
                                         as.numeric(true_hosp$true_value)),

      est_wis = weighted_interval_score(as.numeric(probs0),
                                        spline_quant,
                                        as.numeric(true_hosp$true_value)),
      
      dist_diff = eval_dist(eval_quantile_spline, probs0),
      
      ssp = sum((eval_quantile_spline - probs0)^2),
      
      ssq = sum((quantiles0 - spline_quant)^2),
      
      sap = sum(abs(eval_quantile_spline - probs0)),
      
      saq = sum(abs(quantiles0 - spline_quant)),
      
      n = length(quantiles0)
    )
    
    
   
    
    qkern <- function(p) {qkden(p, quantiles0, kernel = "gaussian")}
    pkern <- function(x) {pkden(x, quantiles0, kernel = "gaussian")}
    rkern <- function(n) {rkden(n, quantiles0, kernel = "gaussian")}
    pukern <- function(x) {pdist(qkern(x))}
    
    eval_quantile_kern <- pkern(quantiles0)
    
    kern_draws <- rkern(sample_size)
    kern_quant <- qkern(probs0)
    
    
    kern_scores0 <- data.frame(
      model = mod,
      date = sub_date, 
      location = loc, 
      horizon = h, 
      method = "kern",
      zeros = TRUE,
      crps = scoringRules::crps_sample(true_hosp$true_value, kern_draws),
      logs = scoringRules::logs_sample(true_hosp$true_value, kern_draws),
      
      true_wis = weighted_interval_score(as.numeric(probs0),
                                         quantiles0,
                                         as.numeric(true_hosp$true_value)),
      
      est_wis = weighted_interval_score(as.numeric(probs0),
                                        kern_quant,
                                        as.numeric(true_hosp$true_value)),
      
      dist_diff = eval_dist(eval_quantile_kern, probs0),
      
      ssp = sum((eval_quantile_kern - probs0)^2),
      
      ssq = sum((quantiles0 - kern_quant)^2),
      
      sap = sum(abs(eval_quantile_kern - probs0)),
      
      saq = sum(abs(quantiles0 - kern_quant)),
      
      n = length(quantiles)
    )
    
    #############################################
    ##############Nothing Removed################
    #############################################
    
    
    qspline <- make_q_fn(probs, quantiles)
    pspline <- make_p_fn(probs, quantiles)
    puspline <- function(x) {pdist(qspline(x))}
    rspline <- make_r_fn(probs, quantiles)
    
    eval_quantile_spline <- pspline(quantiles)
    
    spline_draws <- rspline(sample_size)
    spline_quant <- qspline(probs)
    
    spline_scores <- data.frame(
      model = mod,
      date = sub_date, 
      location = loc, 
      horizon = h, 
      method = "spline",
      zeros = FALSE,
      crps = scoringRules::crps_sample(true_hosp$true_value, spline_draws),
      logs = scoringRules::logs_sample(true_hosp$true_value, spline_draws),
      
      true_wis = weighted_interval_score(as.numeric(probs),
                                          quantiles,
                                          as.numeric(true_hosp$true_value)),
      
      est_wis = weighted_interval_score(as.numeric(probs),
                                        spline_quant,
                                        as.numeric(true_hosp$true_value)),
      
      dist_diff = eval_dist(eval_quantile_spline, probs),
      
      ssp = sum((eval_quantile_spline - probs)^2),
      
      ssq = sum((quantiles - spline_quant)^2),
      
      sap = sum(abs(eval_quantile_spline - probs)),
      
      saq = sum(abs(quantiles - spline_quant)),
      
      n = length(quantiles)
    )
    
    
    qkern <- function(p) {qkden(p, quantiles, kernel = "gaussian")}
    pkern <- function(x) {pkden(x, quantiles, kernel = "gaussian")}
    rkern <- function(n) {rkden(n, quantiles, kernel = "gaussian")}
    pukern <- function(x) {pdist(qkern(x))}
    
    eval_quantile_kern <- pkern(quantiles)
    
    kern_draws <- rkern(sample_size)
    kern_quant <- qkern(probs)
    
    
    kern_scores <- data.frame(
      model = mod,
      date = sub_date, 
      location = loc, 
      horizon = h, 
      method = "kern",
      zeros = FALSE,
      crps = scoringRules::crps_sample(true_hosp$true_value, kern_draws),
      logs = scoringRules::logs_sample(true_hosp$true_value, kern_draws),
      
      true_wis = weighted_interval_score(as.numeric(probs),
                                         quantiles,
                                         as.numeric(true_hosp$true_value)),
      
      est_wis = weighted_interval_score(as.numeric(probs),
                                        kern_quant,
                                        as.numeric(true_hosp$true_value)),
      
      dist_diff = eval_dist(eval_quantile_kern, probs),
      
      ssp = sum((eval_quantile_kern - probs)^2),
      
      ssq = sum((quantiles - kern_quant)^2),
      
      sap = sum(abs(eval_quantile_kern - probs)),
      
      saq = sum(abs(quantiles - kern_quant)),
      
      n = length(quantiles)
    )
    
    # print("Please tell me it gets here!")
    all_scores <- rbind(spline_scores0, kern_scores0,
          spline_scores, kern_scores)
    
    all_scores
    
  }

saveRDS(all_scores, paste0("./flu_forecast_analysis/spl_kern_res/",
                           mod, ".rds"))





