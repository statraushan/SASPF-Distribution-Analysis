############################################################
# Title: Q-Q and P-P Plots for Distribution Comparison
# Dataset: Blood Cancer Data
#
# Author: Shantanu Yadav
# Date: April 2026
############################################################

rm(list = ls(all = TRUE))

############################
# Libraries
############################
library(ggplot2)
library(gridExtra)

############################
# Load and Normalize Data
############################
data <- DataSetsUni::data_bloodcancer
eps <- 1e-6

x <- (data - min(data)) / (max(data) - min(data))
x[x == 0] <- eps
x[x == 1] <- 1 - eps
x <- sort(x)

############################################################
# PDF and CDF Definitions
############################################################

# Custom Distribution
pdf_custom <- function(x, theta) {
  theta * x^(theta - 1) * (x^theta + 1) * exp(x^theta - 1)
}

cdf_custom <- function(x, theta) {
  x^theta * exp(x^theta - 1)
}

# Topp-Leone (TL)
pdf_TL <- function(x, theta) {
  2 * theta * (1 - x) * (x^(theta - 1)) * (2 - x)^(theta - 1)
}

cdf_TL <- function(x, theta) {
  (x * (2 - x))^theta
}

# UNXLD
pdf_UNXLD <- function(x, theta) {
  (theta / 2) * x^(theta - 1) * (1 - theta * log(x))
}

cdf_UNXLD <- function(x, theta) {
  (x^theta) - (theta * (x^theta) * log(x)) / 2
}

# Power Distribution (PD)
pdf_PD <- function(x, theta) {
  theta * x^(theta - 1)
}

cdf_PD <- function(x, theta) {
  x^theta
}

# Kumaraswamy (Kw)
pdf_Kw <- function(x, theta, beta) {
  theta * beta * (x^(theta - 1)) * (1 - x^theta)^(beta - 1)
}

cdf_Kw <- function(x, theta, beta) {
  1 - (1 - x^theta)^beta
}

# Beta Distribution (BD)
pdf_Beta <- function(x, theta, beta) {
  ifelse(
    x > 0 & x < 1,
    (x^(theta - 1)) * ((1 - x)^(beta - 1)) / beta(theta, beta),
    0
  )
}

cdf_Beta <- function(x, theta, beta) {
  pbeta(x, shape1 = theta, shape2 = beta)
}

############################################################
# PDF and CDF Lists with Parameters
############################################################

pdf_list <- list(
  Custom = list(pdf = pdf_custom, theta = 0.597005661),
  TL     = list(pdf = pdf_TL, theta = 1.332726053),
  UNXLD  = list(pdf = pdf_UNXLD, theta = 1.285218764),
  PD     = list(pdf = pdf_PD, theta = 0.909154954),
  Kw     = list(pdf = pdf_Kw, theta = 0.607, beta = 0.578),
  BD     = list(pdf = pdf_Beta, theta = 0.666, beta = 0.59)
)

cdf_list <- list(
  Custom = function(x) cdf_custom(x, 0.597005661),
  TL     = function(x) cdf_TL(x, 1.332726053),
  UNXLD  = function(x) cdf_UNXLD(x, 1.285218764),
  PD     = function(x) cdf_PD(x, 0.909154954),
  Kw     = function(x) cdf_Kw(x, 0.607, 0.578),
  BD     = function(x) cdf_Beta(x, 0.666, 0.59)
)

############################################################
# Empirical CDF for P-P Plots
############################################################
empirical_cdf_vals <- ecdf(x)(x)

############################################################
# Q-Q and P-P Plot Functions
############################################################

qq_plot_func <- function(cdf_func, dist_name) {
  
  u <- ppoints(length(x))
  
  theoretical_q <- sapply(u, function(p) {
    uniroot(function(q) cdf_func(q) - p,
            interval = c(min(x), max(x)))$root
  })
  
  ggplot(
    data.frame(empirical = sort(x), theoretical = sort(theoretical_q)),
    aes(x = theoretical, y = empirical)
  ) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, col = "red") +
    labs(
      title = paste("Q-Q Plot:", dist_name),
      x = "Theoretical Quantiles",
      y = "Empirical Quantiles"
    ) +
    theme_minimal()
}

pp_plot_func <- function(cdf_func, dist_name) {
  
  theoretical_vals <- cdf_func(x)
  
  ggplot(
    data.frame(theoretical = theoretical_vals,
               empirical = empirical_cdf_vals),
    aes(x = theoretical, y = empirical)
  ) +
    geom_point(alpha = 0.7, color = "darkgreen") +
    geom_abline(intercept = 0, slope = 1,
                col = "red", linetype = "dashed") +
    labs(
      title = paste("P-P Plot:", dist_name),
      x = "Theoretical CDF",
      y = "Empirical CDF"
    ) +
    theme_minimal()
}

############################################################
# Generate Q-Q and P-P Plots
############################################################

qq_plots <- lapply(names(cdf_list), function(name) {
  qq_plot_func(cdf_list[[name]], name)
})

pp_plots <- lapply(names(cdf_list), function(name) {
  pp_plot_func(cdf_list[[name]], name)
})

############################################################
# Save Plots
############################################################

pdf("QQ_plot-BC.pdf", width = 6, height = 5)
grid.arrange(grobs = qq_plots, ncol = 2)
dev.off()

pdf("PP_plot-BC.pdf", width = 6, height = 5)
grid.arrange(grobs = pp_plots, ncol = 2)
dev.off()

############################################################
# End of Script
############################################################
