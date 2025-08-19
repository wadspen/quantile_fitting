library(ggplot2)
library(dplyr)


probs <- seq(.025, .975, length.out = 12)
quants <- quantile(rnorm(100), probs = probs)
plot(quants ~ probs)


mod <- lm(quants ~ qnorm(probs))
p <- seq(.02, .98, length.out = 1001)
q <- qnorm(p, coef(mod)[1], coef(mod)[2])
tq <- qnorm(p)


data.frame(probs, quants) %>% 
  ggplot() +
  geom_line(data = data.frame(p, tq), aes(x = p, y = tq), 
            size = 1.2, colour = "lightgrey") +
  geom_line(data = data.frame(p, q), aes(x = p, y = q), 
            size = 1.9, colour = "darkgrey", linetype = "dashed") +
  geom_point(aes(x = probs, y = quants), size = 3.1) +
  xlab("Probability") +
  ylab("Quantile") +
  theme_bw() + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 17))
  
