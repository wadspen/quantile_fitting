#.libPaths("~/rlibs")
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
library(LaplacesDemon)
library(stringr)
library(parallel)
library(doParallel)
library(doMC)
#n.cores <- detectCores()
n.cores <- 62
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

args <- commandArgs()
dist <- args[6]
p <- as.numeric(args[7])
nind <- as.numeric(args[8])
model <- args[9]
num_comps <- as.numeric(args[10])

print(dist); print(p); print(nind); print(model)


mod_loc <- "../stan_models/"

mod_name <- ifelse(model == "clt_shs",
                   "cdf_quantile_normal_mixK_shs.stan",
                   ifelse(model == "clt", "cdf_quantile_normal_mixK.stan",
                          ifelse(model == "ord", "order_normal_mixK_quantiles.stan",
                                 ifelse(model == "clt_hs", 
                                        "cdf_quantile_normal_mixK_hs.stan",
                                        ifelse(model == "ord_hs",
                                               "order_normal_mixK_quantiles_hs.stan",
                                               
                                               ifelse(model == "ord_shs",
                                                      "order_normal_mixK_quantiles_shs.stan",
                                                      ifelse(model == "ind_shs",
                                                             "cdf_ind_quantile_normal_mixK_shs.stan",
                                                             ifelse(model == "clt_sb",
                                                                    "cdf_quantile_normal_mixK_sb.stan",
                                                                    ifelse(model == "ord_sb",
                                                                           "order_normal_mixK_quantiles_sb.stan",
                                                                           ifelse(model == "ind_sb", "cdf_ind_quantile_normal_mixK_sb.stan",
                                                                                  ifelse(model == "ind_hs", "cdf_ind_quantile_normal_mixK_hs.stan",
                                                                                         "cdf_ind_quantile_normal_mixK.stan")))))))))))

mod <- cmdstan_model(stan_file = paste0(mod_loc, mod_name))





burn <- 1000
samples <- 1000
out_s <- 5000
sample_type <- "MCMC"

samp_sizes <- c(50, 150, 500, 1000, 2000, 5000)
levels <- list(
  c(.25, .5, .75),
  c(.025, .25, .5, .95, .975),  
  seq(.1, .9, by = .1),
  c(.05, seq(.1, .9, by = .1), .95),
  c(.025, .05, seq(.1, .9, by = .1), .95, .975),
  c(.01, .025, .05, seq(.1, .9, by = .1), .95, .975, .99),
  seq(.05, .95, by = .05),
  c(.025, seq(.05, .95, by = .05), .975),
  c(.01, .025, seq(.05, .95, by = .05), .975, .99),
  seq(.01, .99, by = .02)	       
)



reps <- 60
n <- samp_sizes[nind]

#p <- 7
#n <- 1000
#replicate <- 4
#p <- 2
#n <- 500
#dist <- "norm"
#replicate <- 2
distance <- foreach(replicate = 1:reps,
                    .packages = c("cmdstanr", "evmix", "distfromq", "EnvStats",
                                  "VGAM", "distr", "dplyr", "LaplacesDemon",
                                  "stringr")
                    ,.errorhandling = "remove"
                    ,.combine = rbind) %dopar% {
                      #foreach(n = samp_sizes, .combine = rbind) %:%
                      #foreach(p = 1:length(levels), .combine = rbind) %dopar% {
                      
                      set.seed(replicate)
                      source("./simulation_functions.R")
                      
                      if (dist == "norm") {
                        samp <- rnorm(n)
                        pdist <- function(x) {pnorm(x)}
                        ddist <- function(x) {dnorm(x)}
                        rdist <- function(n) {rnorm(n)}
                        qdist <- function(p) {qnorm(p)}
                      } else if (dist == "evd") {
                        samp <- revd(n)
                        pdist <- function(x) {pevd(x)}
                        ddist <- function(x) {devd(x)}
                        rdist <- function(n) {revd(n)}
                        qdist <- function(p) {qevd(p)}
                      } else if (dist == "lp") {
                        samp <- rlaplace(n)
                        pdist <- function(x) {plaplace(x)}
                        ddist <- function(x) {dlaplace(x)}
                        rdist <- function(n) {rlaplace(n)}
                        qdist <- function(p) {qlaplace(p)}
                      } else if (dist == "gmix") {
                        pars <- data.frame(mu = c(-1, 1.2),
                                           sigma = c(.9, .6),
                                           weight = c(.35, .65))
                        mdist <- make_gaussian_mixture(pars$mu, pars$sigma, pars$weight)
                        samp <- r(mdist)(n)
                        pdist <- function(x) {p(mdist)(x)}
                        ddist <- function(x) {d(mdist)(x)}
                        rdist <- function(n) {r(mdist)(n)}
                        qdist <- function(p) {q(mdist)(p)}
                      } else if (dist == "tdp") {
                        Ndp <- 5
                        alphadp <- 1
                        
                        v <- rbeta(Ndp - 1, 1, alphadp)
                        pi <- Stick(v)
                        
                        pi <- c(pi, 1 - sum(pi))
                        pi[pi < 0] <- 0
                        while (sum(pi) != 1) {
                          pi <- pi/sum(pi)
                        }
                        
                        mus <- rnorm(Ndp + 1, 0, 2)
                        sigmas <- rep(1, Ndp + 1)
                        sigmas <- rgamma(Ndp + 1, 2)
                        
                        mdist <- make_normmix(mus, sigmas, pi)
                        samp <- r(mdist)(n)
                        pdist <- function(x) {p(mdist)(x)}
                        ddist <- function(x) {d(mdist)(x)}
                        rdist <- function(n) {r(mdist)(n)}
                        qdist <- function(p) {q(mdist)(p)}
                        
                      }
                      
                      probs <- levels[[p]]
                      quantiles <- quantile(samp, probs, type = 2)
                      true_quantiles <- qdist(probs)
                      
                      data <- data.frame(quantile = quantiles, prob = probs,
                                         true_quantile = true_quantiles)
                      stan_data <- make_stan_data(data, size = n, comps = num_comps,
                                                  alph = 1)
                      
                      start <- Sys.time() 
                      fit_mod <- stan_fit_draws(mod, stan_data,
                                                sampler = sample_type, burn = burn, samp = samples,
                                                refresh = 100, out_s = 5000)
                      #fit <- mod$optimize(stan_data, iter = 10000)
                      end <- Sys.time()
                      tdiff <- end - start
                      tdiff <- as.numeric(tdiff, "mins")
                      # fit <- mod$sample(data = stan_data,
                      #                   num_chains = 4,
                      #                   num_cores = 4)
                      
                      
                      
                      fit <- fit_mod[[2]]
                      
                      draws <- fit$draws(format = "df") %>% 
                        filter(if_all(everything(), ~ !is.na(.)))
                      comps <- draws$samp_comp
                      comp_props <- table(comps)/sum(table(comps))
                      comp_props <- comp_props[order(comp_props, decreasing = TRUE)]
                      num_comps95 <- min(which(cumsum(comp_props) >= .95))
                      num_comps98 <- min(which(cumsum(comp_props) >= .98))
                      num_comps99 <- min(which(cumsum(comp_props) >= .99))
                      ncompg1 <- sum(comp_props > .01)
                      
                      # draws <- draws %>%
                      #   # filter(samp_comp %in% as.numeric(names(comp_props))[1:num_comps95])
                      #   filter(samp_comp %in% as.numeric(names(comp_props))[comp_props > .01])
                      
                      
                      drawsQ <- draws %>% 
                        dplyr::select(contains("Q_rep"))
                      
                      
                      data_sum <- data %>% 
                        mutate(low95 = apply(drawsQ, MARGIN = 2, 
                                             FUN = quantile, probs = 0.025),
                               upp95 = apply(drawsQ, MARGIN = 2, 
                                             FUN = quantile, probs = 0.975),
                               mean = apply(drawsQ, MARGIN = 2,
                                            FUN = mean)
                        ) 
                      
                      
                      coverm <- data_sum %>% 
                        mutate(cover = between(true_quantile, 
                                               low95, upp95)) %>% 
                        ungroup() %>% 
                        summarise(m = mean(cover)) %>% 
                        as.numeric()
                      
                      
                      
                      params <- draws %>% 
                        dplyr::select(contains("mus") | contains("sigmas") |
                                        contains("pi"))
                      # m_params <- apply(params, MARGIN = 2, FUN = mean)

		      ds <- floor(round(seq(1, dim(draws)[1], length.out = 200), 0))
                      mix_kl_vec <- c()
                      mix_tv_vec <- c()
                      mix_uwd_vec <- c()
                      for (d in ds) {
                      
                        m_params <- params[d,] #apply(params, MARGIN = 2, FUN = mean)
                        mus <- as.numeric(m_params[str_detect(names(m_params), "mus")])
                        sigmas <- as.numeric(m_params[str_detect(names(m_params), "sigmas")])
                        wts <- as.numeric(m_params[str_detect(names(m_params), "pi")])
                        
                        wts[wts < 0] <- 0
                        wts <- wts/sum(wts)
                        wts <- wts/sum(wts)
                        wts <- wts/sum(wts)
                        wts <- wts/sum(wts)
                        wts <- wts/sum(wts)
                        wts <- wts/sum(wts)
                        
                        
                        sig <- which(wts > 1e-5)
                        
                        wtss <- wts[sig]
                        wtss[which(wtss < 0)] <- 0
                        wtss <- wtss/sum(wtss)
                        wtss <- wtss/sum(wtss)
                        wtss <- wtss/sum(wtss)
                        wtss <- wtss/sum(wtss)
                        wtss <- wtss/sum(wtss)
                        
                        #make_normmix(mus, sigmas, wts)
                        
                        mix_dist <- make_normmix(mus, sigmas, wts)
                        
                        
                        dmix <- function(x) {d(mix_dist)(x)}
                        # qmix <- function(p) {q(mix_dist)(p)}
                        # pmix <- function(x) {p(mix_dist)(x)}
                        # pu <- function(x) {pdist(qmix(x))}
                        # dspline <- make_d_fn(probs, quantiles)
                        # ys <- dspline(x)
                        # ym <- dmix(x)
                        # plot(ym ~ x, type = "l")
                        # lines(yd ~ x, col = "blue")
                        # lines(ys ~ x, col = "red")
                        
                        # d <- d + 1
                        kls <- rdist(samples)
                        py <- ddist(kls)
                        mixx <- dmix(kls)
                        
                        mix_kl_vec[d] <- mean(log(py) - log(mixx))
                        
                        mix_tv_vec[d] <- dens_dist(dmix, ddist)
                        
                        # mix_uwd_vec[d] <- unit_wass_dist(pu, d = 1)
                        # print(d)
                      }
                      
                      mix_kl <- mean(mix_kl_vec, na.rm = TRUE)
                      mix_tv <- mean(mix_tv_vec, na.rm = TRUE)
                      # uwd1 <- mean(mix_uwd_vec, na.rm = TRUE)
                      
                      udraws <- pdist(draws$dist_samp)
                      
                      
                      #make unit ecdfs
                      pu <- function(x) {ecdf(udraws)(x)}
                      
                      
                      
                      uwd1 <- unit_wass_dist(pu, d = 1)
                      
                      
                      qdiff <- sum(abs(quantile(draws$dist_samp, probs = probs) - qdist(probs))^2)
                      
                      data.frame(model = model, n = n, 
                                 numq = length(probs),
                                 dist = dist,
                                 samp_type = sample_type,
                                 rep = replicate,
                                 uwd1 = uwd1, 
                                 kld = mix_kl,
                                 tv = mix_tv,
                                 comp95 = num_comps95, 
                                 comp98 = num_comps98,
                                 comp99 = num_comps99,
                                 compg1 = ncompg1,
                                 coverm = coverm,
                                 qdiff = qdiff,
                                 ftime = tdiff
                      )
                    }


saveRDS(distance, paste0("./simple_res_tm/ncomp", num_comps, "/", 
                         paste(dist, model, p, nind, sep = "_"), ".rds"))
# write.csv(distance, paste0("sim_scores/", dist, "_", sample_type, "_numcomp",
#                            "/", "size", n, "_probs", length(levels[[p]]), "_scores.csv"), row.names = FALSE)



