library(shiny)
library(deSolve)
library(ggplot2)

# parametros do artigo - Tameirao et al. 2022
params <- c(Ka = 0.83, Cl = 0.0039, V1 = 0.31, V2 = 0.13, Q = 0.057)
Tlag <- 7.61

pk_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dAd <- -Ka * Ad
    dAc <- Ka * Ad - (Cl/V1 + Q/V1) * Ac + (Q/V2) * Ap
    dAp <- (Q/V1) * Ac - (Q/V2) * Ap
    list(c(dAd, dAc, dAp), Cp = Ac / V1)
  })
}

ui <- fluidPage(
  titlePanel("PK Florfenicol - Tilapia"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("dose", "Dose (mg/kg):", min = 5, max = 50, value = 15),
      sliderInput("interval", "Intervalo (h):", min = 12, max = 48, value = 24),
      sliderInput("n_doses", "N doses:", min = 1, max = 15, value = 10),
      numericInput("mic", "CIM (ug/mL):", value = 1, min = 0, step = 0.1),
      actionButton("sim", "Simular")
    ),
    mainPanel(
      plotOutput("pk_plot", height = "400px"),
      verbatimTextOutput("info")
    )
  )
)

server <- function(input, output) {
  output$pk_plot <- renderPlot({
    input$sim

    t_end <- (input$n_doses - 1) * input$interval + 72
    times <- seq(0, t_end, by = 0.5)

    events <- data.frame(
      var = "Ad",
      time = seq(0, input$n_doses - 1) * input$interval + Tlag,
      value = input$dose,
      method = "add"
    )

    out <- ode(y = c(Ad = 0, Ac = 0, Ap = 0),
               times = times,
               func = pk_model,
               parms = params,
               events = list(data = events))

    df <- as.data.frame(out)
    df$Cp <- df$Ac / params["V1"]

    p <- ggplot(df, aes(x = time, y = Cp)) +
      geom_line(color = "steelblue", size = 1) +
      geom_hline(yintercept = input$mic, linetype = "dashed", color = "red") +
      labs(x = "Tempo (h)", y = "Cp (ug/mL)",
           title = paste("Dose:", input$dose, "mg/kg")) +
      theme_minimal()

    p
  })

  output$info <- renderPrint({
    input$sim
    cat("App PK basico - florfenicol em tilapias\n")
    cat("Baseado em Tameirao et al. 2022\n")
  })
}

shinyApp(ui, server)
