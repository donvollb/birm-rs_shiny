# Beta Item Response Model With Response Styles (BIRM-RS)
# Shiny app illustrating the BIRM-RS (Vollbracht, Lischetzke, & Henninger, 2026),
# an extension of the Beta Item Response Model (Noel & Dauvier, 2007).
#
# Run locally: shiny::runApp("app.R")
# Required packages: shiny, bslib, shinythemes

library(shiny)
library(bslib)
library(shinythemes)

col_model <- "#cc4778"
col_ref <- "grey55"

# Shape parameters of the BIRM without response styles (reference curve)
birm_shapes <- function(theta, delta, tau) {
  middle <- theta - delta
  c(exp((middle + tau) / 2), exp((-middle + tau) / 2))
}

# ref: optional shape parameters c(m, n) of the same person without response styles
beta_normal <- function(m = 2, n = 2, ref = NULL) {
  y <- seq(0, 1, by = 0.01)
  vals <- dbeta(y, m, n)
  vals_ref <- if (!is.null(ref)) dbeta(y, ref[1], ref[2])
  all_vals <- c(vals, vals_ref)
  plot(
    y,
    vals,
    type = "l",
    col = col_model,
    lwd = 2.5,
    xlim = c(0, 1),
    ylim = c(0, max(all_vals[is.finite(all_vals)]) * if (is.null(ref)) 1 else 1.15),
    xlab = expression(bold(Response ~ Y[ij])),
    ylab = expression(bold(f(Y[ij]))),
    main = bquote(bold(
      "Probability Density Function for B(" *
        .(round(m, 2)) *
        ", " *
        .(round(n, 2)) *
        ")"
    ))
  )
  if (!is.null(ref)) {
    lines(y, vals_ref, col = col_ref, lty = 2, lwd = 2)
    lines(y, vals, col = col_model, lwd = 2.5)
    legend(
      "top",
      horiz = TRUE,
      legend = c("with response styles", "without response styles"),
      col = c(col_model, col_ref),
      lty = c(1, 2),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
}

BIRM <- function(theta, delta, tau, equation = FALSE) {
  middle <- theta - delta
  enumerator_m <- middle + tau
  enumerator_n <- -middle + tau
  m <- exp(enumerator_m / 2)
  n <- exp(enumerator_n / 2)
  if (equation == TRUE) {
    return(round(c(m, n), 2))
  }
  beta_normal(m, n)
}

BIRM_ERS <- function(theta, delta, tau, ers, equation = FALSE, show_ref = FALSE) {
  middle <- theta - delta
  ers_exp <- exp(ers)
  enumerator_m <- ers_exp * middle + tau
  enumerator_n <- -ers_exp * middle + tau
  m <- exp(enumerator_m / 2)
  n <- exp(enumerator_n / 2)
  if (equation == TRUE) {
    return(round(c(m, n), 2))
  }
  ref <- if (show_ref) birm_shapes(theta, delta, tau)
  beta_normal(m, n, ref)
}

BIRM_ARS <- function(theta, delta, tau, ars, x = 1, equation = FALSE, show_ref = FALSE) {
  x <- as.numeric(x)
  middle <- theta - delta + x * ars
  enumerator_m <- middle + tau
  enumerator_n <- -middle + tau
  m <- exp(enumerator_m / 2)
  n <- exp(enumerator_n / 2)
  if (equation == TRUE) {
    return(round(c(m, n), 2))
  }
  ref <- if (show_ref) birm_shapes(theta, delta, tau)
  beta_normal(m, n, ref)
}

BIRM_RS <- function(theta, delta, tau, ers, ars, x = 1, equation = FALSE, show_ref = FALSE) {
  x <- as.numeric(x)
  middle <- theta - delta + x * ars
  ers_exp <- exp(ers)
  enumerator_m <- ers_exp * middle + tau
  enumerator_n <- -ers_exp * middle + tau
  m <- exp(enumerator_m / 2)
  n <- exp(enumerator_n / 2)
  if (equation == TRUE) {
    return(round(c(m, n), 2))
  }
  ref <- if (show_ref) birm_shapes(theta, delta, tau)
  beta_normal(m, n, ref)
}

beta_sample <- function(
  n = 1000,
  delta = 0,
  theta_mu = 0,
  theta_sd = 1,
  tau = 2,
  ars_mu = 0,
  ars_sd = 0.5,
  x = 1,
  ers_mu = 0,
  ers_sd = 0.5,
  show_ref = FALSE
) {
  x <- as.numeric(x)
  theta_vals <- rnorm(n, theta_mu, theta_sd)
  ars_vals <- rnorm(n, ars_mu, ars_sd)
  ers_vals <- rnorm(n, ers_mu, ers_sd)

  # do the above vectorized
  tmp <- sapply(1:n, function(i) {
    BIRM_RS(
      theta = theta_vals[i],
      delta = delta,
      tau = tau,
      ers = ers_vals[i],
      ars = ars_vals[i],
      x = x,
      equation = TRUE
    )
  }) |>
    t() |>
    as.data.frame()

  vec <- rbeta(n, tmp[, 1], tmp[, 2])
  dens <- density(vec)

  # Same persons without response styles
  if (show_ref) {
    shapes_ref <- birm_shapes(theta_vals, delta, tau)
    dens_ref <- density(rbeta(n, shapes_ref[1:n], shapes_ref[(n + 1):(2 * n)]))
  }

  par(mar = c(5.1, 4.5, 4.1, 2.1))

  plot(
    dens,
    col = col_model,
    lwd = 2.5,
    ylim = c(0, max(dens$y, if (show_ref) dens_ref$y) * if (show_ref) 1.15 else 1),
    main = paste0("Density From Simulated Data (n = ", n, ")"),
    xlab = expression(bold(Response ~ Y[ij])),
    ylab = expression(bold(f(hat(Y)[ij])))
  )
  if (show_ref) {
    lines(dens_ref, col = col_ref, lty = 2, lwd = 2)
    lines(dens, col = col_model, lwd = 2.5)
    legend(
      "top",
      horiz = TRUE,
      legend = c("with response styles", "without response styles"),
      col = c(col_model, col_ref),
      lty = c(1, 2),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
}


# Refactor the UI using navbarPage and nav_panel
ui <- navbarPage(
  title = "Beta Item Response Theory With Response Styles: Illustration",
  theme = shinytheme("cerulean"),

  # Beta Normal Tab
  nav_panel(
    "Standard Beta",
    sidebarLayout(
      sidebarPanel(
        withMathJax(HTML("<h4><b>Shape Parameters</b></h4>")),
        sliderInput(
          "m_beta",
          withMathJax(HTML(
            "\\(m_{ij}\\): As \\(m_{ij}\\) increases relative to \\(n_{ij}\\), the distribution shifts to the right."
          )),
          min = 0,
          max = 10,
          value = 2,
          step = 0.5
        ),
        sliderInput(
          "n_beta",
          withMathJax(HTML(
            "\\(n_{ij}\\): As \\(n_{ij}\\) increases relative to \\(m_{ij}\\), the distribution shifts to the left"
          )),
          min = 0,
          max = 10,
          value = 2,
          step = 0.5
        )
      ),
      mainPanel(plotOutput("beta_normal"))
    )
  ),

  # BIRM Tab
  nav_panel(
    "BIRM",
    sidebarLayout(
      sidebarPanel(
        withMathJax(HTML("<h4><b>Person Parameters</b></h4>")),
        sliderInput(
          "theta_birm",
          withMathJax(HTML(
            "\\(\\theta_{i}\\): latent trait (or ability). Represents the trait being measured, as in standard IRT models. Higher values shift the distribution to the right."
          )),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),

        withMathJax(HTML(
          "<hr style='margin: 20px 0; border: 1px solid gray;'>"
        )),

        withMathJax(HTML("<h4><b>Item Parameters</b></h4>")),
        sliderInput(
          "delta_birm",
          withMathJax(HTML(
            "\\(\\delta_{j}\\): item difficulty. Also as in standard IRT models—greater difficulty shifts the distribution to the left."
          )),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),

        sliderInput(
          "tau_birm",
          withMathJax(HTML(
            "\\(\\tau_{j}\\): item dispersion parameter. Conceptual overlap with the discrimination parameter in IRT. Higher values reduce the difference between shape parameters and increase them overall, leading to a narrower distribution slightly shifted toward the center."
          )),
          min = 0,
          max = 3,
          value = 2,
          step = 0.2
        )
      ),
      mainPanel(
        uiOutput("BIRM_equation_title1"),
        uiOutput("BIRM_equation1"),
        uiOutput("BIRM_equation_title2"),
        uiOutput("BIRM_equation2"),
        plotOutput("BIRM")
      )
    )
  ),
  # ERS Only Tab
  nav_panel(
    "ERS only",
    sidebarLayout(
      sidebarPanel(
        HTML("<h4><b>Person Parameters</b></h4>"),
        sliderInput(
          "theta_ers",
          withMathJax(HTML("\\(\\theta_{i}\\): latent trait. Represents the trait being measured, as in standard IRT models. Higher values shift the distribution to the right.")),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "ers_ers",
          withMathJax(HTML("\\(\\eta_{i}^{ERS}\\): ERS (extreme response style) parameter. A positive value amplifies shifts toward the extremes (left or right). A negative value pulls the distribution toward the center, reflecting a middle response style.")),
          min = -1,
          max = 1,
          value = 0,
          step = 0.2
        ),
        withMathJax(HTML(
          "<hr style='margin: 20px 0; border: 1px solid gray;'>"
        )),
        withMathJax(HTML("<h4><b>Item Parameters</b></h4>")),
        sliderInput(
          "delta_ers",
          withMathJax(HTML("\\(\\delta_{j}\\): item difficulty. Also as in standard IRT models—greater difficulty shifts the distribution to the left.")),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "tau_ers",
          withMathJax(HTML("\\(\\tau_{j}\\): item dispersion parameter. Conceptual overlap with the discrimination parameter in IRT. Higher values reduce the difference between shape parameters and increase them overall, leading to a narrower distribution slightly shifted toward the center.")),
          min = 0,
          max = 3,
          value = 2,
          step = 0.2
        ),
        checkboxInput(
          "ref_ers",
          "Show the same person without response styles (grey)",
          TRUE
        )
      ),
      mainPanel(
        uiOutput("ERS_equation_title1"),
        uiOutput("ERS_equation1"),
        uiOutput("ERS_equation_title2"),
        uiOutput("ERS_equation2"),
        plotOutput("ERS")
      )
    )
  ),

  # ARS Only Tab
  nav_panel(
    "ARS Only",
    sidebarLayout(
      sidebarPanel(
        HTML("<h4><b>Person Parameters</b></h4>"),
        sliderInput(
          "theta_ars",
          HTML("\\(\\theta_{i}\\): latent trait. Represents the trait being measured, as in standard IRT models. Higher values shift the distribution to the right."),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "ars_ars",
          HTML("\\(\\eta_{i}^{ARS}\\): ARS (acquiescence response style) parameter. If \\(x_j = 1\\), a positive value shifts the distribution to the right, reflecting agreement. If \\(x_j = -1\\), it shifts to the left due to reverse scoring. A negative value shifts the distribution in the opposite direction, indicating a tendency to disagree."),
          min = -1,
          max = 1,
          value = 0,
          step = 0.2
        ),
        HTML("<hr style='margin: 20px 0; border: 1px solid gray;'>"),
        HTML("<h4><b>Item Parameters</b></h4>"),
        sliderInput(
          "delta_ars",
          HTML("\\(\\delta_{j}\\): item difficulty. Also as in standard IRT models—greater difficulty shifts the distribution to the left."),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "tau_ars",
          HTML("\\(\\tau_{j}\\): item dispersion parameter. Conceptual overlap with the discrimination parameter in IRT. Higher values reduce the difference between shape parameters and increase them overall, leading to a narrower distribution slightly shifted toward the center."),
          min = 0,
          max = 3,
          value = 2,
          step = 0.2
        ),
        radioButtons(
          "x_ars",
          label = HTML("\\(x_{j}\\): Item direction. Coded as 1 for normally scored items and -1 for reverse-scored items."),
          selected = 1,
          choiceNames = c("1: normal-scored", "-1: reverse-scored"),
          choiceValues = c(1, -1)
        ),
        checkboxInput(
          "ref_ars",
          "Show the same person without response styles (grey)",
          TRUE
        )
      ),
      mainPanel(
        uiOutput("ARS_equation_title1"),
        uiOutput("ARS_equation1"),
        uiOutput("ARS_equation_title2"),
        uiOutput("ARS_equation2"),
        plotOutput("ARS")
      )
    )
  ),

  # BIRM-RS Tab
  nav_panel(
    "BIRM-RS",
    sidebarLayout(
      sidebarPanel(
        HTML("<h4><b>Person Parameters</b></h4>"),
        sliderInput(
          "theta_birm_rs",
          HTML("\\(\\theta_{i}\\): latent trait. Represents the trait being measured, as in standard IRT models. Higher values shift the distribution to the right."),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "ers_birm_rs",
          HTML("\\(\\eta_{i}^{ERS}\\): ERS (extreme response style) parameter. A positive value amplifies shifts toward the extremes (left or right). A negative value pulls the distribution toward the center, reflecting a middle response style."),
          min = -1,
          max = 1,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "ars_birm_rs",
          HTML("\\(\\eta_{i}^{ARS}\\): ARS (acquiescence response style) parameter. If \\(x_j = 1\\), a positive value shifts the distribution to the right, reflecting agreement. If \\(x_j = -1\\), it shifts to the left due to reverse scoring. A negative value shifts the distribution in the opposite direction, indicating a tendency to disagree."),
          min = -1,
          max = 1,
          value = 0,
          step = 0.2
        ),
        HTML("<hr style='margin: 20px 0; border: 1px solid gray;'>"),
        HTML("<h4><b>Item Parameters</b></h4>"),
        sliderInput(
          "delta_birm_rs",
          HTML("\\(\\delta_{j}\\): item difficulty. Also as in standard IRT models—greater difficulty shifts the distribution to the left."),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "tau_birm_rs",
          HTML("\\(\\tau_{j}\\): item dispersion parameter. Conceptual overlap with the discrimination parameter in IRT. Higher values reduce the difference between shape parameters and increase them overall, leading to a narrower distribution slightly shifted toward the center."),
          min = 0,
          max = 3,
          value = 2,
          step = 0.2
        ),
        radioButtons(
          "x_birm_rs",
          label = HTML("\\(x_{j}\\): item direction. Coded as 1 for normally scored items and -1 for reverse-scored items."),
          selected = 1,
          choiceNames = c("1: normal-scored", "-1: reverse-scored"),
          choiceValues = c(1, -1)
        ),
        checkboxInput(
          "ref_birm_rs",
          "Show the same person without response styles (grey)",
          TRUE
        )
      ),
      mainPanel(
        uiOutput("BIRM_RS_equation_title1"),
        uiOutput("BIRM_RS_equation1"),
        uiOutput("BIRM_RS_equation_title2"),
        uiOutput("BIRM_RS_equation2"),
        plotOutput("BIRM_RS")
      )
    )
  ),

  # BIRM-RS Simulation Tab
  nav_panel(
    "BIRM-RS Simulation",
    sidebarLayout(
      sidebarPanel(
        withMathJax(HTML("<h4><b>Person Parameters</b></h4>")),
        sliderInput(
          "theta_mu_sim",
          withMathJax(HTML(
            "\\(E(\\theta_{i})\\): expected value of the latent trait"
          )),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "theta_sd_sim",
          withMathJax(HTML(
            "\\(SD(\\theta_{i})\\): standard deviation of the latent trait"
          )),
          min = 0,
          max = 3,
          value = 1,
          step = 0.2
        ),
        withMathJax(HTML(
          "<hr style='margin: 20px 0; border: 1px solid gray;'>"
        )),
        sliderInput(
          "ers_mu_sim",
          withMathJax(HTML(
            "\\(E(\\eta_{i}^{ERS})\\): expected value of the ERS parameter"
          )),
          min = -1,
          max = 1,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "ers_sd_sim",
          withMathJax(HTML(
            "\\(SD(\\eta_{i}^{ERS})\\): standard deviation of the ERS parameter"
          )),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        ),
        withMathJax(HTML(
          "<hr style='margin: 20px 0; border: 1px solid gray;'>"
        )),
        sliderInput(
          "ars_mu_sim",
          withMathJax(HTML(
            "\\(E(\\eta_{i}^{ARS})\\): expected value of the ARS parameter"
          )),
          min = -1,
          max = 1,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "ars_sd_sim",
          withMathJax(HTML(
            "\\(SD(\\eta_{i}^{ARS})\\): standard deviation of the ARS parameter"
          )),
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.1
        ),
        withMathJax(HTML(
          "<hr style='margin: 20px 0; border: 1px solid gray;'>"
        )),
        withMathJax(HTML("<h4><b>Item Parameters</b></h4>")),
        sliderInput(
          "delta_sim",
          withMathJax(HTML("\\(\\delta_{j}\\): item difficulty")),
          min = -3,
          max = 3,
          value = 0,
          step = 0.2
        ),
        sliderInput(
          "tau_sim",
          withMathJax(HTML("\\(\\tau_{j}\\): item dispersion parameter")),
          min = 0,
          max = 3,
          value = 2,
          step = 0.2
        ),
        radioButtons(
          "x_sim",
          label = withMathJax(HTML("\\( x_{j} \\): item direction")),
          selected = 1,
          choiceNames = c("1: normal-scored", "-1: reverse-scored"),
          choiceValues = c(1, -1)
        ),
        checkboxInput(
          "ref_sim",
          "Show the same persons without response styles (grey)",
          TRUE
        )
      ),
      mainPanel(
        uiOutput("simulation_equation_title1"),
        uiOutput("simulation_equation1"),
        uiOutput("simulation_equation_title2"),
        uiOutput("simulation_equation2"),
        uiOutput("simulation_equation_title3"),
        uiOutput("simulation_equation3"),
        plotOutput("simulation")
      )
    )
  )
)

server <- function(input, output) {
  output$beta_normal <- renderPlot({
    beta_normal(m = input$m_beta, n = input$n_beta)
  })

  output$BIRM <- renderPlot({
    BIRM(
      theta = input$theta_birm,
      delta = input$delta_birm,
      tau = input$tau_birm
    )
  })

  # let us do equation1 and equation2 for BIRM and then delete the old equation
  output$BIRM_equation_title1 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations for the Probability Density Function")
    )
  })

  output$BIRM_equation_title2 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations With Values")
    )
  })

  output$BIRM_equation1 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{\\theta_i - \\delta_j + \\tau_j}{2}),$$\n",
        "$$n_{ij} = exp(\\frac{-(\\theta_i - \\delta_j) + \\tau_j}{2})$$\n"
      )
    )
  })

  output$BIRM_equation2 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{",
        input$theta_birm,
        " - ",
        input$delta_birm,
        " + ",
        input$tau_birm,
        "}{2}) = ",
        BIRM(
          theta = input$theta_birm,
          delta = input$delta_birm,
          tau = input$tau_birm,
          equation = TRUE
        )[1],
        ",$$\n",
        "$$n_{ij} = exp(\\frac{-(",
        input$theta_birm,
        " - ",
        input$delta_birm,
        ") + ",
        input$tau_birm,
        "}{2}) = ",
        BIRM(
          theta = input$theta_birm,
          delta = input$delta_birm,
          tau = input$tau_birm,
          equation = TRUE
        )[2],
        "$$\n"
      )
    )
  })

  output$BIRM_RS_equation_title1 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations for the Probability Density Function")
    )
  })

  output$BIRM_RS_equation_title2 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations With Values")
    )
  })

  output$BIRM_RS <- renderPlot({
    BIRM_RS(
      theta = input$theta_birm_rs,
      delta = input$delta_birm_rs,
      tau = input$tau_birm_rs,
      ars = input$ars_birm_rs,
      x = input$x_birm_rs,
      ers = input$ers_birm_rs,
      show_ref = isTRUE(input$ref_birm_rs)
    )
  })

  output$BIRM_RS_equation1 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{exp(\\eta^{ERS}_i) \\cdot (\\theta_i - \\delta_j + x_j \\cdot \\eta^{ARS}_i) + \\tau_j}{2}),$$\n",
        "$$n_{ij} = exp(\\frac{-exp(\\eta^{ERS}_i) \\cdot (\\theta_i - \\delta_j + x_j \\cdot \\eta^{ARS}_i) + \\tau_j}{2})$$\n"
      )
    )
  })

  output$BIRM_RS_equation2 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{exp(",
        input$ers_birm_rs,
        ") \\cdot (",
        input$theta_birm_rs,
        " - ",
        input$delta_birm_rs,
        " + ",
        input$x_birm_rs,
        " \\cdot ",
        input$ars_birm_rs,
        ") + ",
        input$tau_birm_rs,
        "}{2}) = ",
        BIRM_RS(
          theta = input$theta_birm_rs,
          delta = input$delta_birm_rs,
          tau = input$tau_birm_rs,
          ars = input$ars_birm_rs,
          x = input$x_birm_rs,
          ers = input$ers_birm_rs,
          equation = TRUE
        )[1],
        ",$$\n",
        "$$n_{ij} = exp(\\frac{-exp(",
        input$ers_birm_rs,
        ") \\cdot (",
        input$theta_birm_rs,
        " - ",
        input$delta_birm_rs,
        " + ",
        input$x_birm_rs,
        " \\cdot ",
        input$ars_birm_rs,
        ") + ",
        input$tau_birm_rs,
        "}{2}) = ",
        BIRM_RS(
          theta = input$theta_birm_rs,
          delta = input$delta_birm_rs,
          tau = input$tau_birm_rs,
          ars = input$ars_birm_rs,
          x = input$x_birm_rs,
          ers = input$ers_birm_rs,
          equation = TRUE
        )[2],
        "$$\n"
      )
    )
  })

  output$ERS <- renderPlot({
    BIRM_ERS(
      theta = input$theta_ers,
      delta = input$delta_ers,
      tau = input$tau_ers,
      ers = input$ers_ers,
      show_ref = isTRUE(input$ref_ers)
    )
  })

  output$ERS_equation_title1 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations for the Probability Density Function")
    )
  })

  output$ERS_equation_title2 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations With Values")
    )
  })

  # modify the equation to split it (1 and 2)
  output$ERS_equation1 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{exp(\\eta^{ERS}_i) \\cdot (\\theta_i - \\delta_j) + \\tau_j}{2}),$$\n",
        "$$n_{ij} = exp(\\frac{-exp(\\eta^{ERS}_i) \\cdot (\\theta_i - \\delta_j) + \\tau_j}{2})$$\n"
      )
    )
  })

  output$ERS_equation2 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{exp(",
        input$ers_ers,
        ") \\cdot (",
        input$theta_ers,
        " - ",
        input$delta_ers,
        ") + ",
        input$tau_ers,
        "}{2}) = ",
        BIRM_ERS(
          theta = input$theta_ers,
          delta = input$delta_ers,
          tau = input$tau_ers,
          ers = input$ers_ers,
          equation = TRUE
        )[1],
        ",$$\n",
        "$$n_{ij} = exp(\\frac{-exp(",
        input$ers_ers,
        ") \\cdot (",
        input$theta_ers,
        " - ",
        input$delta_ers,
        ") + ",
        input$tau_ers,
        "}{2}) = ",
        BIRM_ERS(
          theta = input$theta_ers,
          delta = input$delta_ers,
          tau = input$tau_ers,
          ers = input$ers_ers,
          equation = TRUE
        )[2],
        "$$\n"
      )
    )
  })

  output$ARS <- renderPlot({
    BIRM_ARS(
      theta = input$theta_ars,
      delta = input$delta_ars,
      tau = input$tau_ars,
      ars = input$ars_ars,
      x = input$x_ars,
      show_ref = isTRUE(input$ref_ars)
    )
  })

  output$ARS_equation_title1 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations for the Probability Density Function")
    )
  })

  output$ARS_equation_title2 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations With Values")
    )
  })

  output$ARS_equation1 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{(\\theta_i - \\delta_j + x_j \\cdot \\eta^{ARS}_i) + \\tau_j}{2}),$$\n",
        "$$n_{ij} = exp(\\frac{-(\\theta_i - \\delta_j + x_j \\cdot \\eta^{ARS}_i) + \\tau_j}{2})$$\n"
      )
    )
  })

  output$ARS_equation2 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{(",
        input$theta_ars,
        " - ",
        input$delta_ars,
        " + ",
        input$x_ars,
        " \\cdot ",
        input$ars_ars,
        ") + ",
        input$tau_ars,
        "}{2}) = ",
        BIRM_ARS(
          theta = input$theta_ars,
          delta = input$delta_ars,
          tau = input$tau_ars,
          ars = input$ars_ars,
          x = input$x_ars,
          equation = TRUE
        )[1],
        ",$$\n",
        "$$n_{ij} = exp(\\frac{-(",
        input$theta_ars,
        " - ",
        input$delta_ars,
        " + ",
        input$x_ars,
        " \\cdot ",
        input$ars_ars,
        ") + ",
        input$tau_ars,
        "}{2}) = ",
        BIRM_ARS(
          theta = input$theta_ars,
          delta = input$delta_ars,
          tau = input$tau_ars,
          ars = input$ars_ars,
          x = input$x_ars,
          equation = TRUE
        )[2],
        "$$\n"
      )
    )
  })

  output$simulation_equation_title1 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations for Data Simulation")
    )
  })

  output$simulation_equation_title2 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Equations With Values")
    )
  })

  output$simulation_equation_title3 <- renderUI({
    tags$div(
      style = "font-weight: bold; font-size: 14px; font-family: sans-serif; text-align: center;",
      tags$b("Parameter Distributions")
    )
  })

  output$simulation_equation1 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{exp(\\eta^{ERS}_i) \\cdot (\\theta_i - \\delta_j + x_j \\cdot \\eta^{ARS}_i) + \\tau_j}{2}),$$\n",
        "$$n_{ij} = exp(\\frac{-exp(\\eta^{ERS}_i) \\cdot (\\theta_i - \\delta_j + x_j \\cdot \\eta^{ARS}_i) + \\tau_j}{2})$$\n"
      )
    )
  })

  output$simulation_equation2 <- renderUI({
    withMathJax(
      paste0(
        "$$m_{ij} = exp(\\frac{exp(\\eta^{ERS}_i) \\cdot (\\theta_i - ",
        input$delta_sim,
        " + ",
        input$x_sim,
        " \\cdot \\eta^{ARS}_i) + ",
        input$tau_sim,
        "}{2}),$$\n",
        "$$n_{ij} = exp(\\frac{-exp(\\eta^{ERS}_i) \\cdot (\\theta_i - ",
        input$delta_sim,
        " + ",
        input$x_sim,
        " \\cdot \\eta^{ARS}_i) + ",
        input$tau_sim,
        "}{2})$$\n"
      )
    )
  })

  output$simulation_equation3 <- renderUI({
    withMathJax(
      paste0(
        "$$\\theta_i \\sim \\mathcal{N}(",
        input$theta_mu_sim,
        ", ",
        input$theta_sd_sim,
        "),$$\n",
        "$$\\eta^{ERS}_i \\sim \\mathcal{N}(",
        input$ers_mu_sim,
        ", ",
        input$ers_sd_sim,
        "),$$\n",
        "$$\\eta^{ARS}_i \\sim \\mathcal{N}(",
        input$ars_mu_sim,
        ", ",
        input$ars_sd_sim,
        ")$$\n"
      )
    )
  })

  output$simulation <- renderPlot({
    beta_sample(
      n = 1000,
      delta = input$delta_sim,
      theta_mu = input$theta_mu_sim,
      theta_sd = input$theta_sd_sim,
      tau = input$tau_sim,
      ars_mu = input$ars_mu_sim,
      ars_sd = input$ars_sd_sim,
      x = input$x_sim,
      ers_mu = input$ers_mu_sim,
      ers_sd = input$ers_sd_sim,
      show_ref = isTRUE(input$ref_sim)
    )
  })
}

shinyApp(ui = ui, server = server)
