
library(ggplot2)
library(tidyr)
library(dplyr)
library(devEMF)
windowsFonts(Calibri = windowsFont("Calibri"))
windowsFonts()

# library(showtext)
# font_add("Calibri", "C:/Windows/Fonts/calibri.ttf")  # Pfad anpassen
# showtext_auto()
# showtext_opts(dpi = 300)

# load data
load("data/df_agg.rda")

sim_plot <- function(df, 
                     items, 
                     type_string = "est_cor", 
                     ylim = NULL,
                     plot_title = "Correlations Between True and Estimated Values", 
                     y_label = "Correlation") {
  
  # 1. Subset & build x
  df_sub <- df %>% 
    filter(item_n1 == items) %>% 
    filter(x_num2 != 0) %>%
    filter(ars_prior == 0.9) %>%
    mutate(
      x = factor(paste0("(", x_num1, "|", x_num2, ")"),
                 levels = unique(paste0("(", x_num1, "|", x_num2, ")")))
    )
  
  
  # 2. Pivot both mean and MCSE into long
  df_long <- df_sub %>%
    pivot_longer(
      cols = matches(paste0("^(mean|MCSE)_", type_string, "_(theta1|theta2|ers|ars)$")),
      names_to  = c(".value", "trait"),
      names_pattern = paste0("^(mean|MCSE)_", type_string, "_(theta1|theta2|ers|ars)$")
    ) %>%
    rename(correlation = mean) %>%
    mutate(
      # clean up trait labels
      trait = recode(trait,
                     theta1 = "Latent Trait 1",
                     theta2 = "Latent Trait 2",
                     ers    = "ERS",
                     ars    = "ARS")
    )

  
  my_cols <- c(
  "Latent Trait 1"  = "#440154",  
  "Latent Trait 2" = "#440154",  
  "ERS"                 = "#35B779",  
  "ARS"                 = "#31688E"  
)
  
  my_ltys <- c(
  "Latent Trait 1"  = "solid",
  "Latent Trait 2" = "dashed",  
  "ARS" = "solid",
  "ERS"                 = "solid"
)
  
  # 3. Plot
p <- ggplot(df_long, aes(x = x, 
                         y = correlation, 
                         colour = trait, 
                         linetype = trait,  
                         group = trait)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = correlation - 1.96 * MCSE,
                    ymax = correlation + 1.96 * MCSE),
                width       = 0.2,
                linetype    = "solid",
                show.legend = FALSE,
              linewidth = 1) +
  labs(
    title  = plot_title,
    x      = "Reverse-Scored Items per Latent Trait",  # einfacher String
    y      = y_label,                                  # dynamisch aus Argument
    colour = "Construct",
    linetype = "Construct"
  ) +
  theme_bw(base_size = 22) +
  theme(
    # text            = element_text(family = "Calibri"),  
    panel.background       = element_rect(fill = "transparent", colour = NA),
    plot.background        = element_rect(fill = "transparent", colour = NA),
    axis.title.x           = element_text(face = "bold"),  # fett X-Achse
    axis.title.y           = element_text(face = "bold"),   # fett Y-Achse
    strip.text          = element_text(face = "bold"),  # fett Facet-Labels
    legend.position = "none",
  ) +
  scale_colour_manual(values = my_cols) +
  scale_linetype_manual(values = my_ltys) +
  facet_grid(n ~ var_ers, labeller = labeller(
    n       = as_labeller(c("600" = "n = 600", "300" = "n = 300")),
    var_ers = as_labeller(c("0.1" = "Var(RS): 0.1",
                            "0.3" = "Var(RS): 0.3",
                            "0.5" = "Var(RS): 0.5"))
  ))

if (!is.null(ylim)) {
  p <- p + ylim(ylim)
}


    return(p)
}

plot_emf <- function(p, filename, width = 17.444, height = 13.902, unit = "cm", ...) {
  ggsave(
    filename = filename,
    plot     = p,
    device   = function(file, ...) devEMF::emf(file = file, family = "Calibri", emfPlusFontToPath = TRUE, ...),
    width    = width, height = height, units = unit,
    dpi      = 300  
  )
}


# Korrelationen

p <- sim_plot(df = df_agg, items = 6, type_string = "est_cor", ylim = c(0.59, 1), plot_title = NULL, y_label = "Correlation")
p
plot_emf(p, "plots/cor_6.emf", width = 18.25, height = 16.91, unit = "cm")

p <- sim_plot(df = df_agg, items = 12, type_string = "est_cor", ylim = c(0.59, 1), plot_title = NULL, y_label = "Correlation")
p
plot_emf(p, "plots/cor_12.emf", width = 18.25, height = 16.91, unit = "cm")

# RMSE

p <- sim_plot(df = df_agg, items = 6, type_string = "est_rmse", ylim = c(0.1, 0.6), plot_title = NULL, y_label = "RMSE")
p
plot_emf(p, "plots/rmse_6.emf", width = 18.25, height = 16.91, unit = "cm")

p <- sim_plot(df = df_agg, items = 12, type_string = "est_rmse", ylim = c(0.1, 0.6), plot_title = NULL, y_label = "RMSE")
p
plot_emf(p, "plots/rmse_12.emf", width = 18.25, height = 16.91, unit = "cm")




