###############################################################################
# Lab 13
# Avery Johnson
###############################################################################

library(tidyverse)
library(xtable)
library(e1071)

###############################################################################
# Question 1
###############################################################################

# part a
zebra <- read.csv("zebrafinches.csv")
view(zebra)

further <- zebra$further
n <- length(further)
skew <- skewness(further)

# t test
t_result <- t.test(further, mu=0, alternative="less")
t_stat <- t_result$statistic

# Gaussian PDF and CDF
fz <- dnorm(t_stat)
Fz <- pnorm(t_stat)

# Edgeworth error approx
error <- (skew / sqrt(n)) * ((2*t_stat^2 + 1) / 6) * fz

(probability <- Fz + error)

# part b
t_vals <- seq(-10, 10, by=0.1)
fz_vals <- dnorm(t_vals)
error_vals <- (skew / sqrt(n)) * ((2*t_vals^2 + 1) / 6) * fz_vals

error_df <- data.frame(t=t_vals, error=error_vals)

ggplot(data = error_df, aes(x=t, y=error)) +
  geom_line(color="red", linewidth=1) +
  xlab("T-Statistic") +
  ylab("Error") +
  ggtitle("First-Order Edgeworth Approximation for Error") +
  theme_bw()

# part c
# critical t value
alpha <- 0.05
error_threshold <- 0.1 * alpha
t_alpha <- qnorm(alpha)
fz_alpha <- dnorm(t_alpha)

# compute required n
numerator <- skew * (2 *t_alpha^2 + 1) * fz_alpha
denominator <- 6 * error_threshold
(n_required <- (numerator / denominator) ^ 2)

###############################################################################
# Question 2
###############################################################################

# part a
closer <- zebra$closer
further <- zebra$further
diff <- zebra$diff

R <- 1000

# closer

sd.closer <- sd(closer)
n.closer <- length(closer)

resamples.closer <- tibble(t_stats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(closer,
                          size = n.closer,
                          replace = T)
  
  resamples.closer$t_stats[i] <- (mean(curr.resample) - 0) / (sd.closer / sqrt(n.closer))
  
}
delta.t.closer <- mean(resamples.closer$t_stats) - 0
resamples.null.closer <- resamples.closer |>
  mutate(t_stats.shifted = t_stats - delta.t.closer)

# mean(resamples.null.closer$t_stats.shifted)

# further

sd.further <- sd(further)
n.further <- length(further)

resamples.further <- tibble(t_stats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(further,
                          size = n.further,
                          replace = T)
  
  resamples.further$t_stats[i] <- (mean(curr.resample) - 0) / (sd.further / sqrt(n.further))
  
}
delta.t.further <- mean(resamples.further$t_stats) - 0
resamples.null.further <- resamples.further |>
  mutate(t_stats.shifted = t_stats - delta.t.further)

# mean(resamples.null.further$t_stats.shifted)

# difference

sd.diff <- sd(diff)
n.diff <- length(diff)

resamples.diff <- tibble(t_stats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(diff,
                          size = n.diff,
                          replace = T)
  
  resamples.diff$t_stats[i] <- (mean(curr.resample) - 0) / (sd.diff / sqrt(n.diff))
  
}
delta.t.diff <- mean(resamples.diff$t_stats) - 0
resamples.null.diff <- resamples.diff |>
  mutate(t_stats.shifted = t_stats - delta.t.diff)

# mean(resamples.null.diff$t_stats.shifted)

# part b

# closer
(p.boot.closer <- mean(resamples.null.closer$t_stats.shifted >= delta.t.closer))
(p.t.closer <- (t.test(x=closer, mu=0, alternative="greater"))$p.value)

# further
(p.boot.further <- mean(resamples.null.further$t_stats.shifted <= delta.t.further))
(p.t.further <- (t.test(x=further, mu=0, alternative="less"))$p.value)

# difference
low <- -delta.t.diff
high <- delta.t.diff

p.low <- mean(resamples.null.diff$t_stats.shifted <= low)
p.high <- mean(resamples.null.diff$t_stats.shifted >= high)
(p.boot.diff <- p.low + p.high)
(p.t.diff <- (t.test(x=diff, mu=0, alternative="two.sided"))$p.value)


# part c
(t_crit.boot.closer <- quantile(resamples.null.closer$t_stats.shifted, 0.05))
(t_crit.t.closer <- qt(0.05, df=length(closer-1)))
  
  
(t_crit.boot.further <- quantile(resamples.null.further$t_stats.shifted, 0.05))
(t_crit.t.further <- qt(0.05, df=length(further-1)))
    
(t_crit.boot.diff <- quantile(resamples.null.diff$t_stats.shifted, 0.05))
(t_crit.t.diff <- qt(0.05, df=length(diff-1)))


# part d
# use resamples
# need this for the flurescence (x bar)

(ci.boot.closer <- quantile(resamples.null.closer$t_stats.shifted, c(0.025, 0.975)))
(ci.low.t.closer <- t.test(x = closer, mu = 0, alternative = "two.sided")$conf.int[1])
(ci.high.t.closer <- t.test(x = closer, mu = 0, alternative = "two.sided")$conf.int[2])


(ci.boot.further <- quantile(resamples.null.further$t_stats.shifted, c(0.025, 0.975)))
(ci.low.t.further <- t.test(x = further, mu = 0, alternative = "two.sided")$conf.int[1])
(ci.high.t.further <- t.test(x = further, mu = 0, alternative = "two.sided")$conf.int[2])


(ci.low.boot.diff <- quantile(resamples.null.diff$t_stats.shifted, c(0.025, 0.975)))
(ci.low.t.diff <- t.test(x = diff, mu = 0, alternative = "two.sided")$conf.int[1])
(ci.high.t.diff <- t.test(x = diff, mu = 0, alternative = "two.sided")$conf.int[2])

###############################################################################
# Question 3
###############################################################################