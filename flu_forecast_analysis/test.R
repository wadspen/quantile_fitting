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
print(locations)
print(length(locations))

