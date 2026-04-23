############################################################
# Title: Simulation Study for Parameter Estimation
# Distribution: Custom Distribution
#
# Author: Shantanu Yadav
# Date: April 2026
############################################################

############################
# Clear Environment
############################
rm(list = ls(all = TRUE))

############################
# Load Required Library
############################
library(lamW)

############################
# Simulation Settings
############################
t <- 1000
n_values <- c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200)

############################
# Initialize Storage Matrices
############################
MLE_n <- matrix(nrow = t, ncol = length(n_values))
colnames(MLE_n) <- paste0("n_", n_values)

MPS_n <- matrix(nrow = t, ncol = length(n_values))
colnames(MPS_n) <- paste0("n_", n_values)

LS_n <- matrix(nrow = t, ncol = length(n_values))
colnames(LS_n) <- paste0("n_", n_values)

CVME_n <- matrix(nrow = t, ncol = length(n_values))
colnames(CVME_n) <- paste0("n_", n_values)

############################
# Random Number Generator
############################
set.seed(133)

generate_random_numbers <- function(theta, n) {
  u <- runif(n)
  w_value <- lambertW0(u * exp(1))
  sqrt_term <- (w_value)^(1 / theta)
  result <- sqrt_term
  return(result)
}

############################
# PDF and CDF Functions
############################
f <- function(x, theta) {
  val <- theta * x^(theta - 1) * (x^theta + 1) * exp(x^theta - 1)
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

F <- function(x, theta) {
  val <- x^theta * exp(x^theta - 1)
  return(pmin(1, pmax(0, val)))
}

############################
# Log-Likelihood Function
############################
LL <- function(par, x) {
  theta <- par
  log_LL <- sum(log(f(x, theta)))
  return(-log_LL)
}

############################
# Maximum Product Spacing
############################
MPS <- function(par, x) {
  theta <- par
  n <- length(x)
  x <- sort(x)
  
  D <- numeric(n + 1)
  D[1] <- F(x[1], theta)
  
  for (j in 2:n) {
    D[j] <- F(x[j], theta) - F(x[j - 1], theta)
  }
  
  D[n + 1] <- 1 - F(x[n], theta)
  
  mps_ll <- (1 / (n + 1)) * sum(log(D))
  return(-mps_ll)
}

############################
# Least Squares Estimator
############################
LS <- function(par, x) {
  theta <- par
  n <- length(x)
  x <- sort(x)
  
  Z <- numeric(n)
  
  for (j in 1:n) {
    Z[j] <- F(x[j], theta) - (j / (n + 1))
  }
  
  ls <- sum(Z^2)
  return(ls)
}

############################
# Cramér-von Mises Estimator
############################
CVME <- function(par, x) {
  theta <- par
  n <- length(x)
  x <- sort(x)
  
  C <- numeric(n)
  
  for (j in 1:n) {
    C[j] <- F(x[j], theta) - ((2 * j - 1) / (2 * n))
  }
  
  cvm <- (1 / (12 * n)) + sum(C^2)
  return(cvm)
}

############################
# Simulation Loop
############################
theta_true <- 0.8   # Alternative values: 2, 4

for (i in 1:length(n_values)) {
  
  n <- n_values[i]
  
  for (k in 1:t) {
    
    gen_data <- generate_random_numbers(theta_true, n)
    x <- gen_data
    
    MLE <- nlm(LL, p = c(theta = 2), x = x, hessian = TRUE)
    MLE_n[k, i] <- MLE$estimate
    
    MPS_est <- nlm(MPS, p = c(theta = 2), x = x, hessian = TRUE)
    MPS_n[k, i] <- MPS_est$estimate
    
    LS_est <- nlm(LS, p = c(theta = 2), x = x, hessian = TRUE)
    LS_n[k, i] <- LS_est$estimate
    
    CVM_est <- nlm(CVME, p = c(theta = 2), x = x, hessian = TRUE)
    CVME_n[k, i] <- CVM_est$estimate
  }
}

############################
# Compute Summary Statistics
############################

# MLE
avg_MLE <- colMeans(MLE_n)
bias_1 <- colMeans(MLE_n - theta_true)
MSE_1 <- colMeans((MLE_n - theta_true)^2)

# MPS
avg_MPS <- colMeans(MPS_n)
bias_2 <- colMeans(MPS_n - theta_true)
MSE_2 <- colMeans((MPS_n - theta_true)^2)

# LS
avg_LSE <- colMeans(LS_n)
bias_3 <- colMeans(LS_n - theta_true)
MSE_3 <- colMeans((LS_n - theta_true)^2)

# CVME
avg_CVME <- colMeans(CVME_n)
bias_5 <- colMeans(CVME_n - theta_true)
MSE_5 <- colMeans((CVME_n - theta_true)^2)

############################
# Store MSE Values
############################
mse_df <- data.frame(
  n = n_values,
  MLE = MSE_1,
  MPS = MSE_2,
  LS = MSE_3,
  CVME = MSE_5
)

############################
# Results Table
############################
results_df <- data.frame(
  n = n_values,
  Avg_MLE = avg_MLE,
  Bias_MLE = bias_1,
  MSE_MLE = MSE_1,
  
  Avg_MPS = avg_MPS,
  Bias_MPS = bias_2,
  MSE_MPS = MSE_2,
  
  Avg_LS = avg_LSE,
  Bias_LS = bias_3,
  MSE_LS = MSE_3,
  
  Avg_CVME = avg_CVME,
  Bias_CVME = bias_5,
  MSE_CVME = MSE_5
)

############################
# Print Results
############################
print(results_df)

############################
# MSE Plot
############################
sample_size <- n_values
mse_ML   <- MSE_1
mse_MPS  <- MSE_2
mse_LSE  <- MSE_3
mse_CVME <- MSE_5

y_max <- max(mse_ML, mse_MPS, mse_LSE, mse_CVME)

pdf("MSE_Plot1.pdf", width = 6, height = 5, pointsize = 12)

plot(
  sample_size, mse_ML,
  type = "b",
  ylim = c(0, y_max),
  xlab = "Sample Size",
  ylab = "MSE",
  col.lab = 4,
  col = "red",
  pch = 19,
  lty = 1,
  axes = FALSE
)

lines(sample_size, mse_MPS, type = "b", col = "green", pch = 17, lty = 1)
lines(sample_size, mse_LSE, type = "b", col = "blue", pch = 15, lty = 1)
lines(sample_size, mse_CVME, type = "b", col = "purple", pch = 18, lty = 1)

axis(1, at = seq(0, 200, by = 20), col.axis = "red")
axis(2, at = seq(0.00, 0.08, by = 0.01), col.axis = "red")

legend(
  100,
  y_max * 0.9,
  legend = c("MSE_ML", "MSE_MPS", "MSE_LSE", "MSE_CVME"),
  col = c("red", "green", "blue", "purple"),
  lty = 1,
  pch = c(19, 17, 15, 18)
)

dev.off()

############################################################
# End of Script
############################################################
