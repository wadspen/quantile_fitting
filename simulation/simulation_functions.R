library(distr)


make_dist_string <- function(par1 = 0, par2 = 1, fam = "Norm") {
  paste0(fam, "(", par1, ",", par2, ")")
}


make_gaussian_mixture <- function(pars) {
  mus <- pars$mu
  sigmas <- pars$sigma
  weights <- pars$weight
  
  dists <- list()
  for (i in 1:nrow(pars)) {
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
  integrate(abs_punit, lower = 0, upper = 1, ppitd_est, d = d,
            subdivisions = 3000)$value^(1/d)
}

make_stan_data <- function(data, size, m = 2, c = 3, sv = 3, nv = 3000,
                           pv = 3) {
  quantiles <- data$quantile
  probs <- data$prob

  dat <- list(
    N = length(quantiles),
    n = size,
    Q = quantiles,
    inv_Phip = qnorm(probs),
    p = probs,
    m = m,
    c = c,
    sv = sv,
    nv = nv,
    pv = pv
  )
  
  dat
}


run_stan_model <- function(mod, data_list, burn = 5000, sample = 5000, 
                           num_chain = 1) {
  
  samps <- mod$sample(data = data_list,
                      iter_warmup = burn,
                      iter_sampling = sample,
                      adapt_delta = .9999,
		      chains = num_chain,
                      refresh = 0)
 
  samps
}


isolate_draws <- function(stan_samps, variable = "dist_samps") {
  draws <- stan_samps$draws(format = "df")
  draws$dist_samps
}

stan_fit_draws <- function(mod, data_list, burn = 5000, sample = 5000,
                      num_chain = 1, variable = "dist_samps") {
	
  samps <- run_stan_model(mod, data_list, burn = 5000, sample = 5000, 
                          num_chain = 1)
  
  draws <- isolate_draws(samps, variable = "dist_samps")
  draws
}


est_quantile <- function(p, draws) {
  quantile(cltnormdraws$dist_samps, p)
}






