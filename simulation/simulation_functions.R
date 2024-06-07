library(distr)


make_dist_string <- function(par1 = 0, par2 = 1, fam = "Norm") {
  paste0(fam, "(", par1, ",", par2, ")")
}


make_gaussian_mixture <- function(mus, sigmas, weights) {
  # mus <- pars$mu
  # sigmas <- pars$sigma
  # weights <- pars$weight
  
  dists <- list()
  for (i in 1:length(mus)) {
    dists[i] <- make_dist_string(mus[i], sigmas[i])
  }
  
  expression <- paste0("UnivarMixingDistribution(", 
        paste(dists, collapse = ","), ",", 
        "mixCoeff=weights", ")")
  
  eval(parse(text = expression))
  
}


rmixdist <- function(n, mixdist) {
  r(mixdist)(n)
}


pmixdist <- function(x, mixdist) {
  p(mixdist)(x)
}


dmixdist <- function(x, mixdist) {
  d(mixdist)(x)
}


qmixdist <- function(p, mixdist) {
  q(mixdist)(x)
}


est_quantiles <- function(sample, probs) {
  quantile(sample, probs)
}


abs_q <- function(p, qdist_est, qdist, d = 1) {
  abs(qdist_est(p) - qdist(p))^d
}

abs_punit <- function(x, ppitd_est, d = 1) {
  abs(ppitd_est(x) - x)^d
}


wass_dist <- function(qdist_est, qdist, d = 1) {
  integrate(abs_q, lower = 0, upper = 1, qdist_est, qdist, d = d,
            subdivisions = 1500)$value^(1/d)
}

unit_wass_dist <- function(ppitd_est, d = 1) {
  (d + 1)*integrate(abs_punit, lower = 0, upper = 1, ppitd_est, d = d,
            subdivisions = 3000)$value
}



psuedo_dist <- function(p1, p2) {
  pi <- min(p1, p2)
  pj <- max(p1, p2)
  num <- pi*(1-pj)
  den <- exp(-.5*(qnorm(pi)^2 + qnorm(pj)^2))
  return(num/den)
}

brownian_corr <- function(p1, p2) {
  pi <- min(p1, p2)
  pj <- max(p1, p2)
  return(pi*(1-pj))
}



make_qcorr <- function(probs) {
  qcorr <- matrix(NA, nrow = length(probs), ncol = length(probs))
  for (i in 1:length(probs)) {
    for (j in 1:length(probs)) {
      # qcorr[i,j] <- psuedo_dist(probs[i], probs[j])
      qcorr[i,j] <- brownian_corr(probs[i], probs[j])
    }
  }
  
  # return(2*pi*qcorr)
  return(qcorr)
  
}



make_stan_data <- function(data, size, comps = 4, m = 2, c = 3, sv = 3,
                           nv = 3000,
                           pv = 3, alpha = 1) {
  quantiles <- data$quantile
  probs <- data$prob

  dat <- list(
    N = length(quantiles),
    n = size,
    Q = quantiles,
    inv_Phip = qnorm(probs),
    QCorr = make_qcorr(probs),
    p = probs,
    n_components = comps,
    m = m,
    c = c,
    sv = sv,
    nv = nv,
    pv = pv, 
    alpha = rep(alpha, comps)
  )
  
  dat
}


run_stan_model <- function(mod, data_list, burn = 5000, sample = 5000, 
                           num_chain = 1, rfresh = 0, sampler = "MCMC",
                           elbo = 1000, grad = 30, out_s = 1000) {
  
  if (sampler == "variational") {
    samps <- mod$variational(data = data_list,
                             elbo_samples = elbo,
                             grad_samples = grad,
                             output_samples = out_s)
  } else {
    samps <- mod$sample(data = data_list,
                        iter_warmup = burn,
                        iter_sampling = sample,
                        adapt_delta = .9999,
	  	      chains = num_chain,
                        refresh = rfresh)
    
  }
 
  samps
}


isolate_draws <- function(stan_samps, variable = "dist_samps") {
  draws <- stan_samps$draws(format = "df")
  draws <- draws$dist_samps
  summary <- stan_samps$summary()
  summary <- summary %>% 
    filter(variable %in% c("mu", "sigma", "n")) %>% 
    dplyr::select(variable, q5, q95) %>% mutate(width = q95 - q5)
  
  return(list(draws = draws, summary = summary))
}

stan_fit_draws <- function(mod, data_list, burn = 5000, sample = 5000,
                      num_chain = 1, variable = "dist_samps", refresh = 0,
                      sampler = "MCMC", elbo = 1000, grad = 30, out_s = 1000) {
	
  samps <- run_stan_model(mod, data_list, burn = burn, sample = sample, 
                          num_chain = 1, rfresh = refresh,
                          sampler = sampler, elbo = elbo, grad = grad,
                          out_s = out_s)
  
  draws <- isolate_draws(samps, variable = "dist_samps")
  draws
}


eval_sum <- function(summary, true_params) {
  sum1 <- summary %>% 
    left_join(true_params, by = "variable") %>% 
    mutate(cover90 = between(truth, q5, q95)) %>% 
    mutate(cover90 = as.numeric(cover90)) %>%
    select(variable, width, cover90) %>% 
    t() %>% 
    row_to_names(1) %>% 
    data.frame() %>%
    mutate_all(as.numeric)
  
  
  
  width <- sum1[1,]
  colnames(width) <- paste("width", colnames(width), sep = "_")
  cover90 <- sum1[2,] 
  colnames(cover90) <- paste("cover90", colnames(cover90), sep = "_")
  
  params_res <- cbind(width, cover90)
  rownames(params_res) <- NULL
  
  if (sum(colnames(params_res) %in% c("width_n", "cover90_n")) == 0) {
    blank_n <- data.frame(width_n = NA, cover90_n = NA)
    params_res <- cbind(params_res, blank_n)
  }
  
  params_res <- params_res %>% 
    select(width_mu, width_sigma, width_n, cover90_mu, cover90_sigma, cover90_n)

  return(params_res)
  
}


est_quantile <- function(p, draws) {
  quantile(cltnormdraws$dist_samps, p)
}






