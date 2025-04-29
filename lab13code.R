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

library(boot)
boot.mean <- function(data, indicies){
  mean(data[indicies])
}

# For closer
boot.closer <- boot(data = closer, statistic = boot.mean, R=10000)
(ci.boot.closer <- boot.ci(boot.closer, type="bca"))
(ci.low.t.closer <- t.test(x = closer, mu = 0, alternative = "two.sided")$conf.int[1])
(ci.high.t.closer <- t.test(x = closer, mu = 0, alternative = "two.sided")$conf.int[2])


boot.further <- boot(data = further, statistic = boot.mean, R=10000)
(ci.boot.further <- boot.ci(boot.further, type="bca"))
(ci.low.t.further <- t.test(x = further, mu = 0, alternative = "two.sided")$conf.int[1])
(ci.high.t.further <- t.test(x = further, mu = 0, alternative = "two.sided")$conf.int[2])


boot.diff <- boot(data = diff, statistic = boot.mean, R=10000)
(ci.boot.diff <- boot.ci(boot.diff, type="bca"))
(ci.low.t.diff <- t.test(x = diff, mu = 0, alternative = "two.sided")$conf.int[1])
(ci.high.t.diff <- t.test(x = diff, mu = 0, alternative = "two.sided")$conf.int[2])


###############################################################################
# Question 3
###############################################################################

# part a

# closer data
R <- 10000
mu0 <- 0
rand.closer <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- closer - mu0
for (i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  rand.closer$xbars[i] <- mean(curr.rand)
}
rand.closer <- rand.closer |>
  mutate(xbars = xbars + mu0) # shifting back

# further data
R <- 10000
mu0 <- 0
rand.further <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- further - mu0
for (i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  rand.further$xbars[i] <- mean(curr.rand)
}
rand.further <- rand.further |>
  mutate(xbars = xbars + mu0) # shifting back

# difference data
R <- 10000
mu0 <- 0
rand.diff <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- diff - mu0
for (i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  rand.diff$xbars[i] <- mean(curr.rand)
}
rand.diff <- rand.diff |>
  mutate(xbars = xbars + mu0) # shifting back

# part b

# closer
delta.closer <- abs(mean(closer) - mu0)
low.closer <- mu0 - delta.closer
high.closer <- mu0 + delta.closer
(p.val.closer <- mean(rand.closer$xbars <= low.closer) +
                mean(rand.closer$xbars >= high.closer))

# further
delta.further <- abs(mean(further) - mu0)
low.further <- mu0 - delta.further
high.further <- mu0 + delta.further
(p.val.further <- mean(rand.further$xbars <= low.further) +
  mean(rand.further$xbars >= high.further))

# difference
delta.diff <- abs(mean(diff) - mu0)
low.diff <- mu0 - delta.diff
high.diff <- mu0 + delta.diff
(p.val.diff <- mean(rand.diff$xbars <= low.diff) +
  mean(rand.diff$xbars >= high.diff))


# part c

# closer
R <- 1000
mu0.iterate <- 0.01
starting.point.closer <- mean(closer)

mu.lower.closer <- starting.point.closer
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.closer <- closer - mu.lower.closer
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.closer *
      sample(x = c(-1, 1),
             size = length(x.shift.closer),
             replace = T)
    
    rand.closer$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.closer <- rand.closer |>
    mutate(xbars = xbars + mu.lower.closer) # shifting back
  
  # p-value 
  delta.closer <- abs(mean(closer) - mu.lower.closer)
  (low.closer <- mu.lower.closer - delta.closer) # mirror
  (high.closer <- mu.lower.closer + delta.closer)   # xbar
  (p.val.closer <- mean(rand.closer$xbars <= low.closer) +
      mean(rand.closer$xbars >= high))
  
  if(p.val.closer < 0.05){
    break
  }else{
    mu.lower.closer <- mu.lower.closer - mu0.iterate
  }
}

mu.upper.closer <- starting.point.closer
repeat{
  rand.closer <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.closer <- closer - mu.upper.closer
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.closer *
      sample(x = c(-1, 1),
             size = length(x.shift.closer),
             replace = T)
    
    rand.closer$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.closer <- rand.closer |>
    mutate(xbars = xbars + mu.upper.closer) # shifting back
  
  # p-value 
  delta.closer <- abs(mean(closer) - mu.upper.closer)
  (low.closer <- mu.upper.closer - delta.closer) # mirror
  (high.closer <- mu.upper.closer + delta.closer)   # xbar
  (p.val.closer <- mean(rand.closer$xbars <= low.closer) +
      mean(rand.closer$xbars >= high))
  
  if(p.val.closer < 0.05){
    break
  }else{
    mu.upper.closer <- mu.upper.closer + mu0.iterate
  }
}

c(mu.lower.closer, mu.upper.closer)

# further
R <- 1000
mu0.iterate <- 0.01
starting.point.further <- mean(further)

mu.lower.further <- starting.point.further
repeat {
  rand.further <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.further <- further - mu.lower.further
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.further *
      sample(x = c(-1, 1),
             size = length(x.shift.further),
             replace = T)
    
    rand.further$xbars[i] <- mean(curr.rand)
  }
  
  rand.further <- rand.further %>%
    mutate(xbars = xbars + mu.lower.further) # shifting back
  
  # p-value 
  delta.further <- abs(mean(further) - mu.lower.further)
  (low.further <- mu.lower.further - delta.further) # mirror
  (high.further <- mu.lower.further + delta.further)   # xbar
  (p.val.further <- mean(rand.further$xbars <= low.further) +
      mean(rand.further$xbars >= high.further))
  
  if(p.val.further < 0.05){
    break
  } else {
    mu.lower.further <- mu.lower.further - mu0.iterate
  }
}

mu.upper.further <- starting.point.further
repeat {
  rand.further <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.further <- further - mu.upper.further
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.further *
      sample(x = c(-1, 1),
             size = length(x.shift.further),
             replace = T)
    
    rand.further$xbars[i] <- mean(curr.rand)
  }
  
  rand.further <- rand.further %>%
    mutate(xbars = xbars + mu.upper.further) # shifting back
  
  # p-value 
  delta.further <- abs(mean(further) - mu.upper.further)
  (low.further <- mu.upper.further - delta.further) # mirror
  (high.further <- mu.upper.further + delta.further)   # xbar
  (p.val.further <- mean(rand.further$xbars <= low.further) +
      mean(rand.further$xbars >= high.further))
  
  if(p.val.further < 0.05){
    break
  } else {
    mu.upper.further <- mu.upper.further + mu0.iterate
  }
}

c(mu.lower.further, mu.upper.further)

# difference
R <- 1000
mu0.iterate <- 0.01
starting.point.diff <- mean(diff)

mu.lower.diff <- starting.point.diff
repeat {
  rand.diff <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.diff <- diff - mu.lower.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.diff *
      sample(x = c(-1, 1),
             size = length(x.shift.diff),
             replace = T)
    
    rand.diff$xbars[i] <- mean(curr.rand)
  }
  
  rand.diff <- rand.diff %>%
    mutate(xbars = xbars + mu.lower.diff) # shifting back
  
  # p-value 
  delta.diff <- abs(mean(diff) - mu.lower.diff)
  (low.diff <- mu.lower.diff - delta.diff) # mirror
  (high.diff <- mu.lower.diff + delta.diff)   # xbar
  (p.val.diff <- mean(rand.diff$xbars <= low.diff) +
      mean(rand.diff$xbars >= high.diff))
  
  if(p.val.diff < 0.05){
    break
  } else {
    mu.lower.diff <- mu.lower.diff - mu0.iterate
  }
}

mu.upper.diff <- starting.point.diff
repeat {
  rand.diff <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.diff <- diff - mu.upper.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.diff *
      sample(x = c(-1, 1),
             size = length(x.shift.diff),
             replace = T)
    
    rand.diff$xbars[i] <- mean(curr.rand)
  }
  
  rand.diff <- rand.diff %>%
    mutate(xbars = xbars + mu.upper.diff) # shifting back
  
  # p-value 
  delta.diff <- abs(mean(diff) - mu.upper.diff)
  (low.diff <- mu.upper.diff - delta.diff) # mirror
  (high.diff <- mu.upper.diff + delta.diff)   # xbar
  (p.val.diff <- mean(rand.diff$xbars <= low.diff) +
      mean(rand.diff$xbars >= high.diff))
  
  if(p.val.diff < 0.05){
    break
  } else {
    mu.upper.diff <- mu.upper.diff + mu0.iterate
  }
}

c(mu.lower.diff, mu.upper.diff)


