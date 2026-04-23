############################################################
# Title: Total Time on Test (TTT) Plots
# Purpose: Examine hazard rate shapes
#
# Author: Shantanu Yadav
# Date: April 2026
############################################################

############################
# Load Required Package
############################
library(DataSetsUni)

############################
# Define epsilon
############################
eps <- 1e-6

############################################################
# Dataset 1: Blood Cancer Data (Increasing Hazard)
############################################################

data <- data_bloodcancer
x <- sort((data - min(data) + eps) / (max(data) - min(data) + 2 * eps))
n <- length(x)

TTT <- c()
for (j in 1:n) {
  s <- 0
  for (i in 1:j) {
    s <- s + x[i]
  }
  TTT[j] <- (s + (n - j) * x[j]) / sum(x)
}

pdf("TTT_Plot1.pdf", width = 8, height = 6, pointsize = 12)

plot(
  0:n / n, c(0, TTT),
  type = "l",
  main = "TTT-Plot",
  col = "red",
  ylab = "T(i/n)",
  xlab = "i/n",
  axes = FALSE
)

abline(a = 0, b = 1)

axis(1, at = seq(0, 1, by = 0.2), col.axis = "red")
axis(2, at = seq(0, 1, by = 0.2), col.axis = "red")

dev.off()

############################################################
# Dataset 2: Increase Rate Data (Increasing Hazard)
############################################################

data <- data_increaserate
x <- sort((data - min(data) + eps) / (max(data) - min(data) + 2 * eps))
n <- length(x)

TTT <- c()
for (j in 1:n) {
  s <- 0
  for (i in 1:j) {
    s <- s + x[i]
  }
  TTT[j] <- (s + (n - j) * x[j]) / sum(x)
}

pdf("TTT_Plot2.pdf", width = 8, height = 6, pointsize = 12)

plot(
  0:n / n, c(0, TTT),
  type = "l",
  main = "TTT-Plot",
  col = "red",
  ylab = "T(i/n)",
  xlab = "i/n",
  axes = FALSE
)

abline(a = 0, b = 1)

axis(1, at = seq(0, 1, by = 0.2), col.axis = "red")
axis(2, at = seq(0, 1, by = 0.2), col.axis = "red")

dev.off()

############################################################
# Dataset 3: Myelogenous Data (Bathtub Hazard)
############################################################

data <- data_Myelogenous
x <- sort((data - min(data) + eps) / (max(data) - min(data) + 2 * eps))
n <- length(x)

TTT <- c()
for (j in 1:n) {
  s <- 0
  for (i in 1:j) {
    s <- s + x[i]
  }
  TTT[j] <- (s + (n - j) * x[j]) / sum(x)
}

pdf("TTT_Plot3.pdf", width = 8, height = 6, pointsize = 12)

plot(
  0:n / n, c(0, TTT),
  type = "l",
  main = "TTT-Plot",
  col = "red",
  ylab = "T(i/n)",
  xlab = "i/n",
  axes = FALSE
)

abline(a = 0, b = 1)

axis(1, at = seq(0, 1, by = 0.2), col.axis = "red")
axis(2, at = seq(0, 1, by = 0.2), col.axis = "red")

dev.off()

############################################################
# End of Script
############################################################
