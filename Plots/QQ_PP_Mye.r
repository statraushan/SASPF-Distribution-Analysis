############################################################
# Title: Q-Q and P-P Plots for Distribution Comparison
# Dataset: Myelogenous Data
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
data <- DataSetsUni::data_Myelogenous
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

# Unit Zeta Distribution (UZD)
pdf_UZD <- function(x, theta) {
  ((theta^3 * x) / ((1 - x)^4 * (theta + 2))) *
    exp(-theta * x / (1 - x))
}

cdf_UZD <- function(x, theta) {
  1 - exp(-theta * x / (1 - x)) *
    (1 + (theta * x * (theta + 2 - 2 * x)) /
       ((theta + 2) * (1 - x)^2))
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

# Unit Transmuted Distribution (UTD)
pdf_UTD <- function(x, theta) {
  theta * ((x^(-theta)) - 1) * x^(-theta - 1) *
    exp(-x^(-theta) + 1)
}

cdf_UTD <- function(x, theta) {
  (x^(-theta)) * exp(-x^(-theta) + 1)
}

############################################################
# PDF and CDF Lists with Parameters
############################################################

pdf_list <- list(
  Custom = list(pdf = pdf_custom, theta = 0.202311323),
  TL     = list(pdf = pdf_TL, theta = 0.383834695),
  UZD    = list(pdf = pdf_UZD, theta = 2.555164215),
  UNXLD  = list(pdf = pdf_UNXLD, theta = 0.462530342),
  PD     = list(pdf = pdf_PD, theta = 0.318597966),
  UTD    = list(pdf = pdf_UTD, theta = 0.401175918)
)

cdf_list <- list(
  Custom = function(x) cdf_custom(x, 0.202311323),
  TL     = function(x) cdf_TL(x, 0.383834695),
  UZD    = function(x) cdf_UZD(x, 2.555164215),
  UNXLD  = function(x) cdf_UNXLD(x, 0.462530342),
  PD     = function(x) cdf_PD(x, 0.318597966),
  UTD    = function(x) cdf_UTD(x, 0.401175918)
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
  theoretical_q <- numeric(length(u))
  
  for (i in seq_along(u)) {
    p <- u[i]
    lower <- min(x)
    upper <- max(x)
    f_lower <- cdf_func(lower) - p
    f_upper <- cdf_func(upper) - p
    
    if (f_lower * f_upper > 0) {
      theoretical_q[i] <- NA
    } else {
      root <- tryCatch(
        uniroot(function(q) cdf_func(q) - p,
                interval = c(lower, upper))$root,
        error = function(e) NA
      )
      theoretical_q[i] <- root
    }
  }
  
  df <- data.frame(empirical = x, theoretical = theoretical_q)
  df <- df[complete.cases(df), ]
  
  ggplot(df, aes(x = theoretical, y = empirical)) +
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

pdf("QQ_plot-M.pdf", width = 6, height = 5)
grid.arrange(grobs = qq_plots, ncol = 2)
dev.off()

pdf("PP_plot-M.pdf", width = 6, height = 5)
grid.arrange(grobs = pp_plots, ncol = 2)
dev.off()

############################################################
# End of Script
############################################################
