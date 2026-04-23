############################################################
# Title: PDF, CDF, Survival and Hazard Function Plots
############################################################

rm(list = ls())

############################
# 1. Generate x values
############################
epsilon <- 1e-6
x <- seq(epsilon, 1 - epsilon, length.out = 1000)

############################
# 2. Define Functions
############################
pdf_fun <- function(x, alpha) {
  alpha * x^(alpha - 1) * (x^alpha + 1) * exp(x^alpha - 1)
}

############################
# 3. Parameter Values
############################
alpha_vals <- c(6, 4, 2, 0.8, 0.6, 0.4)
cols <- c("red", "blue", "green", "purple", "black", "grey")

############################
# 4. Generic Plot Function
############################
plot_distribution <- function(fun, filename, ylab, ylim = NULL,
                              legend_x, legend_y) {
  
  pdf(filename, width = 8, height = 6, pointsize = 12)
  
  plot(
    x, fun(x, alpha_vals[1]),
    type = "l",
    col = cols[1],
    lwd = 1,
    xlab = "x",
    ylab = ylab,
    col.lab = 4,
    axes = FALSE,
    ylim = ylim
  )
  
  for (i in 2:length(alpha_vals)) {
    lines(x, fun(x, alpha_vals[i]), col = cols[i], lwd = 1)
  }
  
  axis(1, at = seq(0, 1, by = 0.2), col.axis = "red")
  axis(2, col.axis = "red")
  
  legend(
    x = legend_x,
    y = legend_y,
    legend = alpha_vals,
    col = cols,
    lty = 1,
    lwd = 1,
    cex = 0.8,
    bty = "n"
  )
  
  dev.off()
}

# PDF Plot
plot_distribution(
  fun = pdf_fun,
  filename = "PDF_Plot.pdf",
  ylab = "Density Function",
  ylim = c(0, 3),
  legend_x = 0.65,
  legend_y = 3
)

############################################################
# End of Script
############################################################
