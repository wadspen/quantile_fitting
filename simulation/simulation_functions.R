library(distr)


eval_dist <- function(ests, probs, p = 1) {
	ests <- c(0, ests, 1)
	probs <- c(0, probs, 1)

	est_diffs <- distancePointToPoint(ests)
	prob_diffs <- distancePointToPoint(probs)

	return((p + 1)*sum(abs(est_diffs - prob_diffs)^p)^(1/p))
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

abs_d <- function(x, ddist_est, ddist, d = 1) {
  abs(ddist_est(x) - ddist(x))^d
}

dens_dist <- function(ddist_est, ddist) {
  (1/2)*integrate(abs_d, lower = -Inf, upper = Inf, ddist_est, ddist, d = 1,
                    subdivisions = 3000)$value
}



norm_dist <- function(p1, p2) {
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

exp_dist <- function(p1, p2) {
  pi <- min(p1, p2)
  return(pi/(1 - pi))
}


exp_dist_test <- function(p1, p2) {
  pi <- min(p1, p2)
  pj <- max(p1, p2)
  return(pi/(1 - pi))
}



make_qcorr <- function(probs, corr = "brown") {
  qcorr <- matrix(NA, nrow = length(probs), ncol = length(probs))
  for (i in 1:length(probs)) {
    for (j in 1:length(probs)) {
      if (corr == "norm") 
        {qcorr[i,j] <- norm_dist(probs[i], probs[j])}
      else if (corr == "brown") {qcorr[i,j] <- brownian_corr(probs[i], probs[j])}
      else if (corr == "exp") {qcorr[i,j] <- exp_dist(probs[i], probs[j])}
    } 
  }
  
  if (corr == "norm") {qcorr <- 2*pi*qcorr}
  return(qcorr)
  
}



make_stan_data <- function(data, size, comps = 4, m = 5, c = 7, sv = 6,
                           nv = 3000,
                           pv = 3000, alpha = 1, cor = "brown") {
  quantiles <- data$quantile
  probs <- data$prob

  dat <- list(
    N = length(quantiles),
    n = size,
    Q = quantiles,
    inv_Phip = qnorm(probs),
    QCorr = make_qcorr(probs, corr = cor),
    p = probs,
    expp = -log(1-probs),
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


isolate_draws <- function(stan_samps, variable = "dist_samp") {
  draws <- stan_samps$draws(format = "df")
  draws <- draws$dist_samp
  summary <- stan_samps$summary(NULL, 
                                posterior::default_summary_measures()[1:4],
                                quantiles = ~ quantile2(., probs = c(0.025, 0.975)),
                                posterior::default_convergence_measures())
  summary <- summary %>% 
    filter(variable %in% c("mu", "sigma", "n")) %>% 
    dplyr::select(variable, q2.5, q97.5) %>% mutate(width = q97.5 - q2.5)
  
  return(list(draws = draws, summary = summary))
}

stan_fit_draws <- function(mod, data_list, burn = 10000, sample = 50000,
                      num_chain = 1, variable = "dist_samp", refresh = 0,
                      sampler = "MCMC", elbo = 1000, grad = 30, out_s = 1000) {
	
  samps <- run_stan_model(mod, data_list, burn = burn, sample = sample, 
                          num_chain = 1, rfresh = refresh,
                          sampler = sampler, elbo = elbo, grad = grad,
                          out_s = out_s)
  
  draws <- isolate_draws(samps, variable = "dist_samp")
  return(list(draws, samps))
}


eval_sum <- function(summary, true_params, quantile_data) {
  sum1 <- summary %>% 
    left_join(true_params, by = "variable") %>% 
    mutate(cover95 = between(truth, q2.5, q97.5)) %>% 
    mutate(cover95 = as.numeric(cover95)) %>%
    mutate(width = q97.5 - q2.5) %>% 
    dplyr::select(variable, width, cover95) %>% 
    t() %>% 
    row_to_names(1) %>% 
    data.frame() %>%
    mutate_all(as.numeric)
  
  
  
  width <- sum1[1,]
  colnames(width) <- paste("width", colnames(width), sep = "_")
  cover95 <- sum1[2,] 
  colnames(cover95) <- paste("cover95", colnames(cover95), sep = "_")
  
  params_res <- cbind(width, cover95)
  rownames(params_res) <- NULL
  
  if (sum(colnames(params_res) %in% c("width_n", "cover95_n")) == 0) {
    blank_n <- data.frame(width_n = NA, cover95_n = NA)
    params_res <- cbind(params_res, blank_n)
  }
  
  params_res <- params_res %>% dplyr::select(contains(true_params[,1])) %>% 
    dplyr::select(-contains("rho"))
    # dplyr::select(width_mu, width_sigma, width_n, cover95_mu, cover95_sigma, cover95_n)
  
  quants_res <- summary %>% 
    filter(str_detect(variable, "Q_rep")) %>% 
    mutate(true_quantile = quantile_data$true_quantile) %>% 
    mutate(width = q97.5 - q2.5, 
           cover = between(true_quantile, q2.5, q97.5)) %>% 
    summarise(width_quant = mean(width),
              cover95_quant =  mean(cover))

  return(cbind(params_res, quants_res))
  
}


est_quantile <- function(p, draws) {
  quantile(cltnormdraws$dist_samp, p)
}






