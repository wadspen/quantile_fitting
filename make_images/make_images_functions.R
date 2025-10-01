library(orderstats)
get_quantile_samps <- function(samps, quantiles, n, n_known = FALSE, n_modeled = FALSE,
                               ind = TRUE, samp_quant = FALSE, order = TRUE,
                               dens_probs = c(.025, .25, .5, .75, .975),
                               true_dist = "lp", model = "cdf") {
  
  
  draws <- samps$draws(format = "df")
  num_samps <- nrow(draws)/10
  draws_s <- draws[1:num_samps*10, ]
  
  all_pis <- draws_s[, str_detect(colnames(draws), "pi")]
  all_mus <- draws_s[, str_detect(colnames(draws), "mus")]
  all_sigmas <- draws_s[, str_detect(colnames(draws), "sigmas")]
  if (ind == TRUE) {all_sigma <- draws_s$sigma}
  if (n_modeled == TRUE) {all_ns <- draws_s$n}
  
  # num_samps <- nrow(all_pis)
  
  
  if (true_dist == "lp") {dx <- seq(-3.5, 3.5, length.out = 1001)} #Laplace bounds
  else if (true_dist == "evd") {dx <- seq(-3, 8, length.out = 1001)} #EVD bounds
  else if (true_dist == "gmix") {dx <- seq(-3, 3, length.out = 1001)}
  m <- 1
  samp_quantiles <- matrix(NA, nrow = num_samps, ncol = length(probs))
  samp_dens <- matrix(NA, nrow = num_samps, ncol = length(dx))
  
  repeat{
    mus <- unlist(all_mus[m,])
    sigmas <- unlist(all_sigmas[m,])
    
    pi <- unlist(all_pis[m,])
    pi[which(pi < 0)] <- 0
    pi <- pi/sum(pi)
    # if (sum(pi) < 1) {
    #   pi[which.max(pi)] <- pi[which.max(pi)] + (1 - sum(pi))
    # } else if (sum(pi) > 1) {
    #   pi[which.min(pi)] <- pi[which.min(pi)] + (1 - sum(pi))
    # }
    if (n_modeled == TRUE) {n <- unlist(all_ns[m])}
    if (ind == TRUE) {sigma <- all_sigma[m]}
    normmix <- UnivarMixingDistribution(Norm(mus[1], sigmas[1]),
                                        Norm(mus[2], sigmas[2]), 
                                        Norm(mus[3], sigmas[3]), 
                                        Norm(mus[4], sigmas[4]),
                                        mixCoeff = pi)
    
    
    
    if (ind == TRUE) {
      fit_samp <- rnorm(length(probs), probs, 1/sigma)
      fit_samp[fit_samp > 1] <- 1
      fit_samp[fit_samp < 0] <- 0
      samp_quantiles[m,] <- q(normmix)(fit_samp)
    }
    if ((n_known == TRUE | n_modeled == TRUE) & samp_quant == FALSE) {
      fit_samp <- mvrnorm(1, probs, (1/n)*stan_data$QCorr)
      fit_samp[fit_samp > 1] <- 1
      fit_samp[fit_samp < 0] <- 0
      samp_quantiles[m,] <- q(normmix)(fit_samp)
    }
    if (samp_quant == TRUE) {
      fit_samp <- r(normmix)(n)
      # fit_samp[fit_samp > 1] <- 1
      # fit_samp[fit_samp < 0] <- 0
      samp_quantiles[m,] <- quantile(fit_samp, probs = probs)
    }
    if (order == TRUE & (n_known == TRUE | n_modeled == TRUE)) {
      osamp <- c()
      for (i in 1:length(probs)){
        osamp[i] <- order_probs(1, probs[i]*n, n)
      }
      samp_quantiles[m,] <- q(normmix)(osamp)
      
    }
    
    samp_dens[m,] <- d(normmix)(dx)
    
    m <- m + 1
    if (m > num_samps) {break}
    # if (m %% 100 == 0) {print(m)}
  }
  
  dens_bounds <- apply(samp_dens, MARGIN = 2, FUN = quantile, 
                       probs = dens_probs)
  dens_bounds <- data.frame(t(dens_bounds))
  colnames(dens_bounds) <- as.character(dens_probs)
  if (true_dist == "lp") {dens_bounds$y <- dlaplace(dx)}
  else if (true_dist == "evd") {dens_bounds$y <- devd(dx)}
  else if (true_dist == "gmix") {dens_bounds$y <- ddist(dx)}
  dens_bounds$x <- dx
  dens_bounds$model <- model
  
  
  quant_bounds <- apply(samp_quantiles, MARGIN = 2, 
                        FUN = quantile, 
                        probs = probs, na.rm = TRUE)
  
  quant_bounds <- data.frame(t(quant_bounds))
  colnames(quant_bounds) <- as.character(probs)
  quant_bounds$prob <- probs
  quant_bounds$quantile <- quantiles
  
  quant_bounds <- quant_bounds %>% 
    mutate(
      cover98 = between(quantile, `0.01`, `0.99`),
      cover95 = between(quantile, `0.025`, `0.975`),
      cover90 = between(quantile, `0.05`, `0.95`),
      cover80 = between(quantile, `0.1`, `0.9`),
      cover70 = between(quantile, `0.15`, `0.85`),
      cover60 = between(quantile, `0.2`, `0.8`),
      cover50 = between(quantile, `0.25`, `0.75`),
      cover40 = between(quantile, `0.3`, `0.7`),
      cover30 = between(quantile, `0.35`, `0.65`),
      cover20 = between(quantile, `0.4`, `0.6`),
      cover10 = between(quantile, `0.45`, `0.55`))
  
  quant_bounds$model <- model
  
  return(list(dens_bounds, quant_bounds))
}


make_dens_plot <- function(data) {
  plot <- data[[1]] %>% 
    ggplot() +
    # geom_histogram(data = data.frame(x = samp), 
    #                aes(x = x, y = ..density..),
    #                fill = "grey") +
    geom_ribbon(aes(x = x, ymin = `0.025`, ymax = `0.975`),
                fill = "pink") + 
    geom_ribbon(aes(x = x, ymin = `0.25`, ymax = `0.75`),
                fill = "red") +
    geom_line(aes(x = x, y = y), size = .8,
              colour = "black") +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=23),
          plot.title=element_text(hjust = .1, margin = margin(t = -4, b=-20)))
  
  return(plot)
}


make_quant_plot <- function(data) {
  plot <- data[[2]] %>% 
    # filter(between(prob, .025, .975)) %>%
    ggplot() +
    geom_ribbon(aes(x = prob, ymin = `0.025`, ymax = `0.975`),
                fill = "pink") + 
    geom_ribbon(aes(x = prob, ymin = `0.25`, ymax = `0.75`),
                fill = "red") +
    # geom_line(aes(x = prob, y = `0.5`), 
    #           colour = "purple", size = 1) +
    geom_point(aes(x = prob, y = quantile),
               colour = "black", size = 1.3) +
    # coord_cartesian(ylim=c(-2.6, 5.5)) + #EVD
    # coord_cartesian(ylim=c(-5, 5)) + #Laplace
    coord_cartesian(ylim=c(-2.5, 3)) + #gmix
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=23),
          plot.title=element_text(hjust = .1, margin = margin(t = 0, b=-20)))
  
  return(plot)
  
}
