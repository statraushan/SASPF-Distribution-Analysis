############################################################
# Title: Distribution Fitting and Goodness-of-Fit Analysis
# Dataset: Myelogenous Data
############################################################

############################
# Clear Environment
############################
rm(list = ls(all = TRUE))

############################
# Load Libraries
############################
library(numDeriv)
library(goftest)
library(DataSetsUni)

############################
# Data Preprocessing
############################
eps <- 1e-6
data <- data_Myelogenous

x <- (data - min(data) + eps) / (max(data) - min(data) + 2 * eps)
n <- length(x)

############################################################
# PDF & CDF Definitions
############################################################

# Custom Distribution
pdf_custom <- function(x, theta) {
  theta * x^(theta - 1) * (1 + x^theta) * exp(x^theta - 1)
}

cdf_custom <- function(x, theta) {
  x^theta * exp(x^theta - 1)
}

# Topp-Leone Distribution
pdf_TL <- function(x, theta) {
  2 * theta * (1 - x) * x^(theta - 1) * (2 - x)^(theta - 1)
}

cdf_TL <- function(x, theta) {
  (x * (2 - x))^theta
}

# Unit Zeghdoudi Distribution
pdf_UZD <- function(x, theta) {
  val <- ((theta^3 * x) / ((1 - x)^4 * (theta + 2))) *
    exp(-theta * x / (1 - x))
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_UZD <- function(x, theta) {
  val <- (1 - exp(-theta * x / (1 - x))) *
    (1 + (theta * x * (theta + 2 - 2 * x)) /
       ((theta + 2) * (1 - x)^2))
  return(pmin(1, pmax(0, val)))
}

# Unit New XLindley Distribution
pdf_UNXLD <- function(x, theta) {
  (theta / 2) * x^(theta - 1) * (1 - theta * log(x))
}

cdf_UNXLD <- function(x, theta) {
  x^theta - (theta * x^theta * log(x)) / 2
}

# Power Distribution
pdf_PD <- function(x, theta) {
  theta * x^(theta - 1)
}

cdf_PD <- function(x, theta) {
  x^theta
}

# Unit Teissier Distribution
pdf_UTD <- function(x, theta) {
  val <- theta * ((x^(-theta)) - 1) * x^(-theta - 1) *
    exp(-x^(-theta) + 1)
  val[!is.finite(val) | val <= 0] <- .Machine$double.eps
  return(val)
}

cdf_UTD <- function(x, theta) {
  val <- x^(-theta) * exp(-x^(-theta) + 1)
  return(pmin(1, pmax(0, val)))
}

############################################################
# Distribution List
############################################################
dist_list <- list(
  Custom = list(pdf = pdf_custom, cdf = cdf_custom),
  TL     = list(pdf = pdf_TL, cdf = cdf_TL),
  UZD    = list(pdf = pdf_UZD, cdf = cdf_UZD),
  UNXLD  = list(pdf = pdf_UNXLD, cdf = cdf_UNXLD),
  PD     = list(pdf = pdf_PD, cdf = cdf_PD),
  UTD    = list(pdf = pdf_UTD, cdf = cdf_UTD)
)

############################################################
# Negative Log-Likelihood
############################################################
negloglik <- function(theta, pdf, x) {
  
  if (theta <= 0) return(Inf)
  
  val <- pdf(x, theta)
  
  if (any(!is.finite(val)) || any(val <= 0)) return(Inf)
  
  return(-sum(log(val)))
}

############################################################
# Fitting Function
############################################################
fit_model <- function(pdf, cdf, x) {
  
  opt <- tryCatch({
    optim(
      par = 2,
      fn = negloglik,
      pdf = pdf,
      x = x,
      method = "L-BFGS-B",
      lower = 1e-4
    )
  }, error = function(e) NULL)
  
  if (is.null(opt) || opt$convergence != 0) return(NULL)
  
  theta_hat <- opt$par
  LL <- -opt$value
  
  ############################
  # Hessian for Standard Error
  ############################
  H <- tryCatch({
    numDeriv::hessian(function(t) negloglik(t, pdf, x), theta_hat)
  }, error = function(e) NA)
  
  SE <- if (!any(is.na(H))) sqrt(1 / H) else NA
  
  ############################
  # Goodness-of-Fit Test
  ############################
  cdf_fn <- function(q) pmin(1, pmax(0, cdf(q, theta_hat)))
  
  ks <- tryCatch(ks.test(x, cdf_fn), error = function(e) NULL)
  
  ############################
  # Information Criteria
  ############################
  k <- 1
  AIC  <- -2 * LL + 2 * k
  AICc <- AIC + (2 * k^2 + 2 * k) / (length(x) - k - 1)
  BIC  <- -2 * LL + log(length(x)) * k
  
  return(data.frame(
    Theta = theta_hat,
    SE = SE,
    LL = LL,
    AIC = AIC,
    AICc = AICc,
    BIC = BIC,
    KS_stat = if (!is.null(ks)) ks$statistic else NA,
    KS_p    = if (!is.null(ks)) ks$p.value else NA
  ))
}

############################################################
# Apply to All Models
############################################################
results <- lapply(names(dist_list), function(name) {
  
  cat("Fitting:", name, "\n")
  
  res <- fit_model(
    pdf = dist_list[[name]]$pdf,
    cdf = dist_list[[name]]$cdf,
    x = x
  )
  
  if (is.null(res)) {
    return(data.frame(
      Distribution = name,
      Theta = NA, SE = NA, LL = NA,
      AIC = NA, AICc = NA, BIC = NA,
      KS_stat = NA, KS_p = NA
    ))
  }
  
  cbind(Distribution = name, res)
})

############################################################
# Final Table
############################################################
final_results <- do.call(rbind, results)

numeric_cols <- sapply(final_results, is.numeric)
final_results[numeric_cols] <- round(final_results[numeric_cols], 4)

print(final_results)

############################################################
# End of Script
############################################################
