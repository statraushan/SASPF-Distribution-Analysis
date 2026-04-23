############################################################
# Title: Distribution Fitting and Goodness-of-Fit Analysis
# Dataset: Increase Rate Data (Second data: Blood Cancer)
############################################################

############################
# Clear Environment
############################
rm(list = ls(all = TRUE))

############################
# Load Libraries
############################
library(MASS)
library(goftest)
library(DataSetsUni)
library(ADGofTest)

############################
# Load and Normalize Data
############################
eps <- 1e-6
data <- data_bloodcancer
#data_increaserate ## second data: data_bloodcancer
x <- (data - min(data) + eps) / (max(data) - min(data) + 2 * eps)
n <- length(x)

############################################################
# PDF and CDF Definitions
############################################################

# Custom Distribution
pdf_custom <- function(x, theta) {
  val <- theta * x^(theta - 1) * (x^theta + 1) * exp(x^theta - 1)
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_custom <- function(x, theta) {
  val <- x^theta * exp(x^theta - 1)
  return(pmin(1, pmax(0, val)))
}

# Topp-Leone Distribution
pdf_TL <- function(x, theta) {
  val <- 2 * theta * (1 - x) * x^(theta - 1) * (2 - x)^(theta - 1)
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_TL <- function(x, theta) {
  val <- (x * (2 - x))^theta
  return(pmin(1, pmax(0, val)))
}

# Unit New XLindley Distribution
pdf_UNXLD <- function(x, theta) {
  val <- (theta / 2) * x^(theta - 1) * (1 - theta * log(x))
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_UNXLD <- function(x, theta) {
  val <- x^theta - (theta * x^theta * log(x)) / 2
  return(pmin(1, pmax(0, val)))
}

# Power Distribution
pdf_PD <- function(x, theta) {
  val <- theta * x^(theta - 1)
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_PD <- function(x, theta) {
  val <- x^theta
  return(pmin(1, pmax(0, val)))
}

# Kumaraswamy Distribution
pdf_Kw <- function(x, theta, beta) {
  val <- theta * beta * x^(theta - 1) * (1 - x^theta)^(beta - 1)
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_Kw <- function(x, theta, beta) {
  val <- 1 - (1 - x^theta)^beta
  return(pmin(1, pmax(0, val)))
}

############################################################
# Distribution Lists
############################################################

single_param_distributions <- list(
  Custom = list(pdf = pdf_custom, cdf = cdf_custom),
  ToppLeone = list(pdf = pdf_TL, cdf = cdf_TL),
  UnitNewXLindley = list(pdf = pdf_UNXLD, cdf = cdf_UNXLD),
  PowerDist = list(pdf = pdf_PD, cdf = cdf_PD)
)

############################################################
# Fit Single-Parameter Distributions
############################################################

results_single <- lapply(names(single_param_distributions), function(name) {
  
  dist <- single_param_distributions[[name]]
  
  fit <- tryCatch({
    fitdistr(x, dist$pdf, start = list(theta = 2), lower = 0.001)
  }, error = function(e) NULL)
  
  if (is.null(fit)) {
    return(data.frame(
      Distribution = name, Theta = NA, SE = NA, LL = NA,
      AIC = NA, AICc = NA, BIC = NA, KS_stat = NA, KS_p = NA
    ))
  }
  
  theta_hat <- fit$estimate
  se_hat <- fit$sd
  LL <- fit$loglik
  k <- length(theta_hat)
  
  AIC  <- -2 * LL + 2 * k
  AICc <- AIC + (2 * k^2 + 2 * k) / (n - k - 1)
  BIC  <- -2 * LL + log(n) * k
  
  cdf_fn <- function(q) dist$cdf(q, theta_hat)
  ks <- tryCatch(ks.test(x, cdf_fn), error = function(e) NULL)
  
  data.frame(
    Distribution = name,
    Theta = round(theta_hat, 4),
    SE = round(se_hat, 4),
    LL = round(LL, 4),
    AIC = round(AIC, 4),
    AICc = round(AICc, 4),
    BIC = round(BIC, 4),
    KS_stat = if (!is.null(ks)) round(ks$statistic, 4) else NA,
    KS_p = if (!is.null(ks)) round(ks$p.value, 4) else NA
  )
})

############################################################
# Fit Multi-Parameter Distributions
############################################################

multi_param_results <- list()

# Kumaraswamy
fit_Kw <- tryCatch({
  fitdistr(x, pdf_Kw, start = list(theta = 2, beta = 2),
           lower = c(0.001, 0.001))
}, error = function(e) NULL)

if (!is.null(fit_Kw)) {
  theta_hat <- fit_Kw$estimate
  se_hat <- fit_Kw$sd
  LL <- fit_Kw$loglik
  k <- length(theta_hat)
  AIC <- -2 * LL + 2 * k
  AICc <- AIC + (2 * k^2 + 2 * k) / (n - k - 1)
  BIC <- -2 * LL + log(n) * k
  
  cdf_fn <- function(q) cdf_Kw(q, theta_hat["theta"], theta_hat["beta"])
  ks <- tryCatch(ks.test(x, cdf_fn), error = function(e) NULL)
  
  multi_param_results[["Kumaraswamy"]] <- data.frame(
    Distribution = "Kumaraswamy",
    Theta = paste0(round(theta_hat, 4), collapse = ", "),
    SE = paste0(round(se_hat, 4), collapse = ", "),
    LL = round(LL, 4),
    AIC = round(AIC, 4),
    AICc = round(AICc, 4),
    BIC = round(BIC, 4),
    KS_stat = if (!is.null(ks)) round(ks$statistic, 4) else NA,
    KS_p = if (!is.null(ks)) round(ks$p.value, 4) else NA
  )
}

# Beta Distribution
fit_Beta <- tryCatch({
  fitdistr(x, dbeta, start = list(shape1 = 2, shape2 = 2),
           lower = c(0.001, 0.001))
}, error = function(e) NULL)

if (!is.null(fit_Beta)) {
  theta_hat <- fit_Beta$estimate
  se_hat <- fit_Beta$sd
  LL <- fit_Beta$loglik
  k <- length(theta_hat)
  AIC <- -2 * LL + 2 * k
  AICc <- AIC + (2 * k^2 + 2 * k) / (n - k - 1)
  BIC <- -2 * LL + log(n) * k
  
  cdf_fn <- function(q) pbeta(q, shape1 = theta_hat["shape1"], shape2 = theta_hat["shape2"])
  ks <- tryCatch(ks.test(x, cdf_fn), error = function(e) NULL)
  
  multi_param_results[["Beta"]] <- data.frame(
    Distribution = "Beta",
    Theta = paste0(round(theta_hat, 4), collapse = ", "),
    SE = paste0(round(se_hat, 4), collapse = ", "),
    LL = round(LL, 4),
    AIC = round(AIC, 4),
    AICc = round(AICc, 4),
    BIC = round(BIC, 4),
    KS_stat = if (!is.null(ks)) round(ks$statistic, 4) else NA,
    KS_p = if (!is.null(ks)) round(ks$p.value, 4) else NA
  )
}

############################################################
# Combine and Print Results
############################################################

final_results <- do.call(rbind, c(results_single, multi_param_results))
print(format(final_results, nsmall = 4))

############################################################
# End of Script
############################################################
