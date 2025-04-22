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

probability <- Fz + error

# part b
t_vals <- seq(-10, 10, by=0.1)
fz_vals <- dnorm(t_vals)
error_vals <- (skew / sqrt(n)) * ((2*t_vals^2 + 1) / 6) * fz_vals

error_df <- data.frame(t=t_vals, error=error_vals)

ggplot(error_df, aes(x=t, y=error)) +
  geom_line(color="red", linewidth=1) +
  xlab("t-statistic") +
  ylab("error") +
  ggtitle("First-Order Edgeworth Approximation Error") +
  theme_bw()

# part c


###############################################################################
# Question 2
###############################################################################





###############################################################################
# Question 3
###############################################################################