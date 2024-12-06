library(cmdstanr)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(stringr)
library(distr)
library(VGAM)
library(distr)

mod <- cmdstan_model(stan_file = 
                       './stan_models/fit_gaussian_mix5_quantiles.stan')

indmod <- cmdstan_model(stan_file = 
                          './stan_models/fit_gaussian_mix5_ind_quantiles.stan')

comp <- 5
n <- 400
probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

yln <- rlnorm(n)
ylp <- rlaplace(n)
yun <- runif(n)

quantilesln <- quantile(yln, probs = probs)
quantileslp <- quantile(ylp, probs = probs)
quantilesun <- quantile(yun, probs = probs)

datln <- list(
  N = length(quantilesln),
  Q = quantilesln,
  p = probs,
  n_components = 5,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

datlp <- list(
  N = length(quantileslp),
  Q = quantileslp,
  p = probs,
  n_components = 5,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

datun <- list(
  N = length(quantilesun),
  Q = quantilesun,
  p = probs,
  n_components = 5,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

sampsln <- mod$sample(data = datln,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)

indsampsln <- indmod$sample(data = datln,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)

drawsindln <- indsampsln$draws(format = "df")

drawsln <- sampsln$draws(format = "df")
diagnosticsln <- sampsln$diagnostic_summary()
summaryln <- sampsln$summary()
listln <- list(draws = drawsln, diagnostics = diagnosticsln, 
               summary = summaryln)

saveRDS(listln, "test_ln.rds")

sampslp <- mod$sample(data = datlp,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)

drawslp <- sampslp$draws(format = "df")
diagnosticslp <- sampslp$diagnostic_summary()
summarylp <- sampslp$summary()

saveRDS(listlp, "test_lp.rds")

sampsun <- mod$sample(data = datun,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)

drawsun <- sampsun$draws(format = "df")
diagnosticsun <- sampsun$diagnostic_summary()
summaryun <- sampsun$summary()

saveRDS(listun, "test_un.rds")


##########################################################################
############Two component Gaussian mixture fit with same model############
##########################################################################


mod2 <- cmdstan_model(stan_file = 
                       './stan_models/fit_gaussian_mix2_quantiles.stan')

comp <- 2
mdist <- UnivarMixingDistribution(Norm(1, 1.1), 
                                  Norm(10.5, 1.1),
                                  mixCoeff = c(.4, .6))

ymd <- r(mdist)(n)

quantilesmd <- quantile(ymd, probs = probs)
datmd <- list(
  N = length(quantilesmd),
  Q = quantilesmd,
  p = probs,
  n_components = comp,
  m = 2,
  c = 3,
  sv = 3,
  tv = 1,
  nv = 3000
)

sampsmd <- mod2$sample(data = datmd,
                      iter_warmup = 5000,
                      iter_sampling = 2000,
                      chains = 1)

drawsmd <- sampsmd$draws(format = "df")
diagnosticsmd <- sampsmd$diagnostic_summary()
summarymd <- sampsmd$summary()
listmd <- list(draws = drawsmd, diagnostics = diagnosticsmd, 
               summary = summarymd)

saveRDS(listln, "test_md.rds")


















