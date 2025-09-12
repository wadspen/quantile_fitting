source("../simulation/simulation_functions.R")
library(dplyr)
library(scoringRules)
library(stringr)
library(evalcast)
library(lubridate)
library(evmix)
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)


args <- commandArgs()
mod <- args[6]
print(mod)


forc_loc <- "../../../FluSight-forecast-hub/model-output/"
hosp_loc <- paste0(forc_loc, "../target-data/target-hospital-admissions.csv")
mod_loc <- "../model-fits-ind/"
models <- list.files(forc_loc)
sub_dates <- substr(list.files(paste0(forc_loc, 
                                "FluSight-baseline/")),
                    1, 10)
horizons <- 0:3
get_loc_file <- list.files(paste0(forc_loc, "FluSight-baseline/"))[4]
get_loc_forc <- read.csv(paste0(forc_loc, "FluSight-baseline/", get_loc_file))
locations <- unique(get_loc_forc$location)

hosp_data <- read.csv(hosp_loc) %>% 
  #select(-X) %>% 
  filter(year(date) >= 2023)
#
#mod <- "FluSight-lop_norm"
#h <- 2
#loc <- "17"
#sub_date <- "2023-10-28"

 all_forecasts <- forecasts <- foreach(sub_date = sub_dates,
                                       .packages = c("distr", "dplyr", 
                                                     "stringr", "scoringRules",
                                                     "evalcast", "LaplacesDemon")
                                       ,.errorhandling = "remove"
                                       ,.combine = rbind) %:%
   #foreach(sub_date = sub_dates, .combine = rbind) %:%
     foreach(loc = locations) %:%
       foreach(h = horizons, .combine = rbind) %dopar% {
       #start1 <- Sys.time()
       if (dir.exists(paste0(mod_loc, mod, "/coverage")) == FALSE) {
          dir.create(paste0(mod_loc, mod, "/coverage"))
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
        
        
        num_samps <- 650 
	draws_s <- draws[sample(nrow(draws), num_samps, replace = TRUE), ]
         
        all_pis <- draws_s[, str_detect(colnames(draws), "pi")]
        all_mus <- draws_s[, str_detect(colnames(draws), "mus")]
        all_sigmas <- draws_s[, str_detect(colnames(draws), "sigmas")]
        all_ns <- draws_s$n
        
        
        #start <- Sys.time() 
        m <- 1
        samp_quantiles <- matrix(NA, nrow = num_samps, ncol = length(probs))
        repeat{
          num <- sample(num_samps, 1)
          mus <- unlist(all_mus[num,])
          sigmas <- unlist(all_sigmas[num,])
          pis <- unlist(all_pis[num,])
          pis[which(pis < 0)] <- 0
          if (sum(pis) < 1) {
            pis[which.max(pis)] <- pis[which.max(pis)] + (1 - sum(pis))
          } else if (sum(pis) > 1) {
            pis[which.min(pis)] <- pis[which.min(pis)] + (1 - sum(pis))
          }
          n <- unlist(all_ns[num])
          n <- ceiling(n) 
          normmix <- UnivarMixingDistribution(Norm(mus[1], sigmas[1]),
                                              Norm(mus[2], sigmas[2]), 
                                              Norm(mus[3], sigmas[3]), 
                                              Norm(mus[4], sigmas[4]), 
                                              mixCoeff = pis)
          samp <- r(normmix)(n)

	  #samp <- rnormm(n, pis, mus, sigmas)
          samp_quantiles[m,] <- quantile(samp, probs = probs)
          
          m <- m + 1; print(m)
          if (m > num_samps) {break}
        }
        #end <- Sys.time()
        
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
                 cover10 = between(quantile, `0.45`, `0.55`),
	  	 wid98 = `0.99` - `0.01`,
	  	 wid95 = `0.975` - `0.025`,
	  	 wid90 = `0.95` - `0.05`,
	  	 wid80 = `0.9` - `0.1`,
	  	 wid70 = `0.85` - `0.15`,
	  	 wid60 = `0.8` - `0.2`,
	  	 wid50 = `0.75` - `0.25`,
	  	 wid40 = `0.7` - `0.3`,
	  	 wid30 = `0.65` - `0.35`,
	  	 wid20 = `0.6` - `0.4`,
	  	 wid10 = `0.55` - `0.45`)
	#end1 <- Sys.time()
		saveRDS(quant_bounds, paste0(mod_loc, mod,
				 "/coverage/", sub_date, "-", loc,
				 "-", h, "-", mod, ".rds"))
##
 #       loop_time <- end - start
 #       all_time <- end1 - start1
 #       quant_bounds$loop_time <- loop_time
 #       quant_bounds$all_time <- all_time
 #       saveRDS(quant_bounds, "test_cover.rds")
    
     
       }
























