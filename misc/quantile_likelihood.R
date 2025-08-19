library(likelihoodExplore)
library(distr)
normal_censor_likelihood <- function(par, fx, quantiles, probs, n) {
  mu <- par[1]
  sigma <- par[2]
  pdiff <- probs[2:length(probs)] - probs[1:(length(probs) - 1)]
  
  left <- (pnorm(quantiles[1], mu, sigma))^(min(probs)*n)
  right <- (1 - pnorm(quantiles[length(quantiles)], mu, sigma))^
    ((1 - max(probs))*n)
  pqdiff <- (pnorm(quantiles[2:length(quantiles)], mu, sigma) - 
    pnorm(quantiles[1:(length(quantiles) - 1)], mu, sigma))
  
  for (i in 1:length(pdiff)) {
    pqdiff[i] <- pqdiff[i]^(pdiff[i]*n)
  }
  
  # log(left) + log(right) + sum(log(pqdiff))# - sum(log(fx))
  (left)*(right)*
    prod((pqdiff))  #- prod(fx)
}
dist <- UnivarMixingDistribution(Norm(-1,.9),
                                 Norm(3, 1.1),
                                 mixCoeff = c(.35, .65))

probs <- seq(.00000001, .999999, length.out = 1000)
probs <- c(.01, .025, seq(.05, .95, by = .05), .975, .99)
probs <- seq(.5, .95, by = .05)
probs <- seq(.01, .99, by = .01)
# prob <- c(.025, seq(.05,.95, by = .05), .975)
samp <- rnorm(400, 15, 3)
quantiles <- quantile(samp, probs)
breaks <- quantile(samp, probs)
# quantiles <- q(dist)(probs)

get_quantile_derivative <- function(quantiles, probs) {
  pdiff <- probs[2:length(probs)] - probs[1:(length(probs) - 1)]
  qdiff <- quantiles[2:length(quantiles)] - 
    quantiles[1:(length(quantiles) - 1)]
  
  qdiff/pdiff
}


fx <- 1/get_quantile_derivative(quantiles, probs)
optim(par = c(10, 2), normal_censor_likelihood, quantiles = quantiles, fx = fx)$par
sd(quantiles)
prod(fx*qdiff)
x <- pnorm(probs, 15, 3)
plot(fx ~ quantiles[-length(quantiles)], type = "l")
hist(samp, breaks = 23)


# prob <- c(.025, seq(.05,.95, by = .05), .975)
probs <- round(seq(.01, .99, by = .01), 2)
# probs <- round(seq(.05, .95, by = .05), 2)
# probs <- seq(0, 1, length.out = 21)
probs <- c(.1,seq(.25, .75, by = .05),.9)
# probs <- c(seq(.25, .75, by = .05))
M <- 1000
m <- 1
mmus <- c()
msigmas <- c()
# repeat{
samp <- rnorm(1000, 3, 2)
quantiles <- quantile(samp, probs)

# samp <- rnorm(21, 15, 3)
# quantiles <- samp[order(samp)]

mut <- 12
sigmat <- 5
mus <- c()
sigmas <- c()
for (i in 1:10000) {
  if(is.na(mut)) {print("prob with mut"); break}
  mustar <- rnorm(1,0,2) + mut
  if(is.na(mustar)) {print("prob with mustar"); break}
  aprob <- min((
                normal_censor_likelihood(c(mustar, sigmat), 
                           fx, quantiles, probs, n) * 
                              dnorm(mustar, 0, 10) *
                              dnorm(mut, mustar, 1)) /
                 
              (
                normal_censor_likelihood(c(mut, sigmat), 
                              fx, quantiles, probs, n) * 
                              dnorm(mut, 0, 10) *
                              dnorm(mustar, mut, 1)), 
                1)
  
  if(is.na(aprob)) {print("prob with aprob"); break}
  accept <- rbinom(1,1,aprob)
  mut <- ifelse(accept, mustar, mut)
  
  ##########################################
  ##########################################
  sigmastar <- rgamma(1, 2, .2)
  aprob <- min((
    normal_censor_likelihood(c(mut, sigmastar), 
                             fx, quantiles, probs, n) * 
      dgamma(sigmastar, 2, .5) *
      dgamma(sigmat, 1, 4)) /
      
      (
        normal_censor_likelihood(c(mut, sigmat), 
                                 fx, quantiles, probs, n) * 
          dgamma(sigmat, 2, .5) *
          dgamma(sigmastar, 1, 4)), 
    1)
  
  if(is.na(aprob)) {print("prob with aprob"); break}
  accept <- rbinom(1,1,aprob)
  sigmat <- ifelse(accept, sigmastar, sigmat)
  
  
  
  
  # print(accept)
  mus[i] <- mut
  sigmas[i] <- sigmat
}
mean(mus[5000:10000]); mean(sigmas[5000:10000]); sd(quantiles)
# mmus[m] <- mean(mus[2500:5000])
# msigmas[m] <- mean(sigmas[2500:5000])
# print(m)
# m <- m + 1
# if (m > M) {break}
# }


sampx <- c()
for (i in 1:2500) {
  mu_samp <- sample(mus[2500:5000], 1)
  sigma_samp <- sample(sigmas[2500:5000], 1)
  
  sampx[i] <- rnorm(1, mu_samp, sigma_samp)
  
}

















