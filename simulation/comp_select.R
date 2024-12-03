library(mixR)
library(distr)
library(VGAM)
library(EnvStats)

n <- 5000
comps <- 2:6
aicomps <- c()
bicomps <- c()
for (i in 1:100) {
	aics <- c()
	bics <- c()
	samp <- rlaplace(n)
	for (comp in 1:length(comps)) {
		mod <- mixfit(samp, ncomp = comps[comp])
		aics[comp] <- mod$aic
		bics[comp] <- mod$bic
	}
	aicomps[i] <- comps[which.min(aics)]
	bicomps[i] <- comps[which.min(bics)]
	print(i)
}


saveRDS(data.frame(aicomp = aicomps, bicomp = bicomps), "lp_comps.rds")

comps <- 2:6
aicomps <- c()
bicomps <- c()
for (i in 1:100) {
	aics <- c()
	bics <- c()
	samp <- revd(n)
	for (comp in 1:length(comps)) {
		mod <- mixfit(samp, ncomp = comps[comp])
		aics[comp] <- mod$aic
		bics[comp] <- mod$bic
	}
	aicomps[i] <- comps[which.min(aics)]
	bicomps[i] <- comps[which.min(bics)]
	print(i)
}


saveRDS(data.frame(aicomp = aicomps, bicomp = bicomps), "evd_comps.rds")

pars <- data.frame(mu = c(-1, 1.2),
                         sigma = c(.9, .6),
                         weight = c(.35, .65))

mix <- UnivarMixingDistribution(Norm(-1, .9),
				Norm(1.2, .6),
				mixCoeff = c(.35, .65))


comps <- 2:6
aicomps <- c()
bicomps <- c()
for (i in 1:100) {
	aics <- c()
	bics <- c()
	samp <- r(mix)(n)
	for (comp in 1:length(comps)) {
		mod <- mixfit(samp, ncomp = comps[comp])
		aics[comp] <- mod$aic
		bics[comp] <- mod$bic
	}
	aicomps[i] <- comps[which.min(aics)]
	bicomps[i] <- comps[which.min(bics)]
}


saveRDS(data.frame(aicomp = aicomps, bicomp = bicomps), "gmix_comps.rds")
