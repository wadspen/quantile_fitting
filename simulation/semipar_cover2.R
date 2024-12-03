source("./simulation_functions.R")
#setwd("./simulation/")
library(cmdstanr)
library(distr)
library(distfromq)
library(evmix)
library(dplyr)
library(VGAM)
library(EnvStats)
library(tidyr)
library(orderstats)
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
dist <- args[6]
p <- as.numeric(args[7])
nind <- as.numeric(args[8])

print(dist); print(p); print(nind)


samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
levels <- list(
  seq(.1, .9, by = .1),
  c(.05, seq(.1, .9, by = .1), .95),
  c(.025, .05, seq(.1, .9, by = .1), .95, .975),
  c(.01, .025, .05, seq(.1, .9, by = .1), .95, .975, .99),
  seq(.05, .95, by = .05),
  c(.025, seq(.05, .95, by = .05), .975),
  c(.01, .025, seq(.05, .95, by = .05), .975, .99),
  seq(.01, .99, by = .02)	       
		)

models <- c("cltn", "ordn", "clt", "ord", "ind")


reps <- 1000
n <- samp_sizes[nind]


#p <- 2
#n <- 500
#dist <- "lp"
#replicate <- 2
coverage <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                  "VGAM", "distr", "dplyr", "orderstats")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
 #foreach(n = samp_sizes, .combine = rbind) %:%
 #foreach(p = 1:length(levels), .combine = rbind) %dopar% {
    
    source("./coverage_functions.R") 
    source("./simulation_functions.R")

    probs <- levels[[p]]
    if (dist == "evd") {
      quantiles <- qevd(probs)
    } else if (dist == "lp") {
      samp <- rlaplace(n)
      quantiles <- qlaplace(probs)
    } else if (dist == "gmix") {
      pars <- data.frame(mu = c(-1, 1.2),
                         sigma = c(.9, .6),
                         weight = c(.35, .65))
      mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
      quantiles <- q(mdist)(probs)
    }
    
       
    data <- data.frame(quantile = quantiles, prob = probs)
    stan_data <- make_stan_data(data, size = n, comps = 4)
    
    print("one")    
            
    file <- paste0("sim_draws/", dist, "_", models[1], "_rep", replicate, 
			"_size", n, "_probs", length(probs), ".rds") 
    file_rds <- readRDS(file)
    quant_bounds1 <- get_fit_coverage(file_rds[[1]], true_quantiles = quantiles,
			             true_probs = probs, n = n, n_known = TRUE,
				     ind = FALSE, QCorr = stan_data$QCorr)
    print("two")
    file <- paste0("sim_draws/", dist, "_", models[2], "_rep", replicate, 
			"_size", n, "_probs", length(probs), ".rds") 
    file_rds <- readRDS(file)
    quant_bounds2 <- get_fit_coverage(file_rds[[1]], true_quantiles = quantiles,
			             true_probs = probs, n = n, n_known = TRUE,
				     order = TRUE, ind = FALSE, QCorr = stan_data$QCorr) 
   

    print("three")
    file <- paste0("sim_draws/", dist, "_", models[3], "_rep", replicate, 
			"_size", n, "_probs", length(probs), ".rds") 
    file_rds <- readRDS(file)
    quant_bounds3 <- get_fit_coverage(file_rds[[1]], true_quantiles = quantiles,
			             true_probs = probs, n = n, n_modeled = TRUE,
				     ind = FALSE,
				     QCorr = stan_data$QCorr) 
 
    print("four")
    file <- paste0("sim_draws/", dist, "_", models[4], "_rep", replicate, 
			"_size", n, "_probs", length(probs), ".rds") 
    file_rds <- readRDS(file)
    quant_bounds4 <- get_fit_coverage(file_rds[[1]], true_quantiles = quantiles,
			             true_probs = probs, n = n, n_modeled = TRUE,
				     order = TRUE, ind = FALSE,
				     QCorr = stan_data$QCorr)


    print("five")
    file <- paste0("sim_draws/", dist, "_", models[5], "_rep", replicate, 
			"_size", n, "_probs", length(probs), ".rds") 
    file_rds <- readRDS(file)
    quant_bounds5 <- get_fit_coverage(file_rds[[1]], true_quantiles = quantiles,
			             true_probs = probs, n = n, ind = TRUE,
				     QCorr = stan_data$QCorr) 


    print("one more")
#    file <- paste0("sim_draws/", dist, "_", mod[6], "_rep", replicate, 
#			"_size", n, "_probs", length(probs), ".rds") 
#    file_rds <- readRDS(file)
#    quant_bounds6 <- get_fit_coverage(file_rds[[1]], true_quantiles = quantiles,
#			             true_probs = probs, n = n, QCorr = stan_data$QCorr)  



    quant_bounds <- rbind(quant_bounds1, quant_bounds2, 
			  quant_bounds3, quant_bounds4,
			  quant_bounds5)

    quant_bounds$model <- rep(models, each = length(probs))
    quant_bounds$rep <- replicate
    quant_bounds
}


probs <- levels[[p]]
cover_path <- paste0("sim_coverage/", dist,  
			"/size", n, "_probs", length(probs), ".rds") 

saveRDS(coverage, cover_path)


#write.csv(distance, paste0("sim_scores/", dist, "/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)


