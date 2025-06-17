
library(ggplot2)
library(tidyr)
library(dplyr)
library(devEMF)
windowsFonts(CambriaMath = windowsFont("Cambria Math"))
windowsFonts()

beta_normal <- function(m = 2, n = 2) {
  vals <- dbeta(seq(0, 1, by = 0.01), m, n)
  plot(seq(0, 1, by = 0.01), vals, type = "l", xlim = c(0, 1), 
       xlab = expression(bold(Y[ij])), 
       ylab = expression(bold(f(Y[ij]))), 
       main = bquote(bold("Probability Density Function for B(" * .(round(m, 2)) * ", " * .(round(n, 2)) * ")")))
}


BIRM_RS <- function(theta = 0, delta = 0, tau = 2, ers = 0, ars = 0, x = 1, equation = TRUE) {
  x <- as.numeric(x)
  middle <- theta - delta + x * ars
  ers_exp <- exp(ers)
  enumerator_m <- ers_exp * middle + tau
  enumerator_n <- -ers_exp * middle + tau
  m <- exp(enumerator_m / 2)
  n <- exp(enumerator_n / 2)
  if (equation == TRUE) return(round(c(m, n), 2))
  beta_normal(m, n)
}

# put the part beginning with df in a function

# plot_demo <- function(vals1, vals2, vals3) {
#   df <- data.frame(x = seq(0, 1, by = 0.01), 
#                    y1 = dbeta(seq(0, 1, by = 0.01), vals1[1], vals1[2]),
#                    y2 = dbeta(seq(0, 1, by = 0.01), vals2[1], vals2[2]),
#                    y3 = dbeta(seq(0, 1, by = 0.01), vals3[1], vals3[2]))
  
#   plot(x = df$x, y = df$y1, type = "l", col = "black", 
#        xlim = c(0, 1), ylim = c(0, 2.5), 
#        xlab = expression(bold(Response ~ Y[ij])), 
#        ylab = expression(bold(Probability ~ Density ~ f(Y[ij]))) 
#        # main = "Probability Density Functions"
#       )
  
#   lines(x = df$x, y = df$y2, col = "#f89540", lwd = 2)
#   lines(x = df$x, y = df$y3, col = "#cc4778", lwd = 2)
  
# }


plot_demo <- function(vals1, vals2, vals3, leg_labels = c("0", "Positive", "Negative")) {
  # 1. build the same data.frame
  x <- seq(0, 1, by = 0.01)
  df <- data.frame(
    x  = x,
    y1 = dbeta(x, vals1[1], vals1[2]),
    y2 = dbeta(x, vals2[1], vals2[2]),
    y3 = dbeta(x, vals3[1], vals3[2])
  )
  
  # 2. pivot to long format
  df_long <- pivot_longer(df,
                          cols      = starts_with("y"),
                          names_to  = "dist",
                          values_to = "density")
  
  # 3. make the plot
  ggplot(df_long, aes(x = x, y = density, colour = dist)) +
    geom_line(linewidth = 2) +
    
    # 4. manually set your colours & line-widths
    scale_colour_manual(
      values = c(y1 = "black", y2 = "#f89540", y3 = "#cc4778")
      # labels = c(y1 = leg_labels[1], y2 = leg_labels[2], y3 = leg_labels[3])
    ) +
    
    # 5. axes, limits, and math expressions
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2.5)) +
    labs(
      x     = expression(bold(Response ~ Y[ij])),
      # y     = expression(bold(Probability ~ Density ~ f(Y[ij]))),
      y     = expression(bold(Probability ~ Density)),
      colour = ""
    ) +
    theme_bw(base_size = 23) + 
    theme(legend.position = "none",
      panel.background  = element_rect(fill = "transparent", colour = NA),
      plot.background   = element_rect(fill = "transparent", colour = NA)) 
}

plot_emf <- function(p, filename, width = 12.46, height = 9.93, unit = "cm") {
  ggsave(
    filename = filename,
    plot     = p,
    device   = function(file, ...) devEMF::emf(file = file, family = "Calibri", emfPlusFontToPath = TRUE, ...),
    width    = width, height = height, units = unit,
    dpi      = 300  
  )
}



# THETA 
vals1 <- BIRM_RS(theta = 0)
vals2 <- BIRM_RS(theta = 1)
vals3 <- BIRM_RS(theta = -1)
p <- plot_demo(vals1, vals2, vals3)
p
plot_emf(p, "plots/theta.emf")

# DELTA
vals1 <- BIRM_RS(delta = 0)
vals2 <- BIRM_RS(delta = 1)
vals3 <- BIRM_RS(delta = -1)
plot_demo(vals1, vals2, vals3) -> p
p
plot_emf(p, "plots/delta.emf")

# ERS
vals1 <- BIRM_RS(theta = 0.5, ers = 0)
vals2 <- BIRM_RS(theta = 0.5, ers = 0.7)
vals3 <- BIRM_RS(theta = 0.5, ers = -0.7)

plot_demo(vals1, vals2, vals3) -> p
p
plot_emf(p, "plots/ers.emf")

# ARS and X (identical)
vals1 <- BIRM_RS(ars = 0)
vals2 <- BIRM_RS(ars = 0.5)
vals3 <- BIRM_RS(ars = -0.5)

plot_demo(vals1, vals2, vals3) -> p
p
plot_emf(p, "plots/ars_x.emf")


# TAU
vals1 <- BIRM_RS(tau = 2)
vals2 <- BIRM_RS(tau = 2.5)
vals3 <- BIRM_RS(tau = 1.5)

plot_demo(vals1, vals2, vals3) -> p
p
plot_emf(p, "plots/tau.emf")



