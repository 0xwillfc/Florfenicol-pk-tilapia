# =============================================================================
# Modelo Farmacocinetico de Florfenicol em Tilapias (Oreochromis niloticus)
# Baseado em: Tameirao et al. (2022) - Rev Inv Vet Peru 33(6): e22433
# Modelo: oral, absorcao de 1a ordem, Tlag, 2 compartimentos, eliminacao linear
# Nesta versao, a temperatura fica fixa na referencia publicada (28 C).
# =============================================================================

library(shiny)
library(deSolve)
library(ggplot2)
library(dplyr)
library(bslib)

# =============================================================================
# CONSTANTES E PARAMETROS
# =============================================================================

PK_PARAMS <- list(
  Tlag_pop = 7.61,
  Ka_pop = 0.83,
  Cl_pop = 0.0039,
  V1_pop = 0.31,
  V2_pop = 0.13,
  Q_pop = 0.057,
  T_ref = 28
)

PK_LIMITS <- list(
  Tlag_min = 0, Tlag_max = 48,
  Ka_min = 0.001, Ka_max = 50,
  Cl_min = 0.0001, Cl_max = 10,
  V1_min = 0.001, V1_max = 10,
  V2_min = 0.001, V2_max = 10,
  Q_min = 0, Q_max = 10
)

MAIN_CURVE_COLOR <- "#1F4E79"
DOSE_COLORS <- c(
  "10 mg/kg" = "#0072B2",
  "15 mg/kg" = "#D55E00",
  "20 mg/kg" = "#009E73"
)

# =============================================================================
# FUNCOES AUXILIARES
# =============================================================================

validate_pk_params <- function(pk_params, lim = PK_LIMITS) {
  validate(
    need(!is.na(pk_params$Tlag_pop) &&
         pk_params$Tlag_pop >= lim$Tlag_min &&
         pk_params$Tlag_pop <= lim$Tlag_max,
         paste0("Tlag deve estar entre ", lim$Tlag_min, " e ", lim$Tlag_max, " h.")),
    need(!is.na(pk_params$Ka_pop) &&
         pk_params$Ka_pop >= lim$Ka_min &&
         pk_params$Ka_pop <= lim$Ka_max,
         paste0("Ka deve estar entre ", lim$Ka_min, " e ", lim$Ka_max, " 1/h.")),
    need(!is.na(pk_params$Cl_pop) &&
         pk_params$Cl_pop >= lim$Cl_min &&
         pk_params$Cl_pop <= lim$Cl_max,
         paste0("Cl deve estar entre ", lim$Cl_min, " e ", lim$Cl_max, " L/h/kg.")),
    need(!is.na(pk_params$V1_pop) &&
         pk_params$V1_pop >= lim$V1_min &&
         pk_params$V1_pop <= lim$V1_max,
         paste0("V1 deve estar entre ", lim$V1_min, " e ", lim$V1_max, " L/kg.")),
    need(!is.na(pk_params$V2_pop) &&
         pk_params$V2_pop >= lim$V2_min &&
         pk_params$V2_pop <= lim$V2_max,
         paste0("V2 deve estar entre ", lim$V2_min, " e ", lim$V2_max, " L/kg.")),
    need(!is.na(pk_params$Q_pop) &&
         pk_params$Q_pop >= lim$Q_min &&
         pk_params$Q_pop <= lim$Q_max,
         paste0("Q deve estar entre ", lim$Q_min, " e ", lim$Q_max, " L/h/kg."))
  )
}

time_above_mic <- function(time, Cp, mic) {
  n <- length(time)
  if (n < 2) return(0)

  above <- Cp >= mic
  t_total <- 0

  for (i in seq_len(n - 1)) {
    t1 <- time[i]
    t2 <- time[i + 1]
    c1 <- Cp[i]
    c2 <- Cp[i + 1]
    a1 <- above[i]
    a2 <- above[i + 1]

    if (a1 && a2) {
      t_total <- t_total + (t2 - t1)
    } else if (a1 && !a2 && c2 != c1) {
      t_cross <- t1 + (mic - c1) / (c2 - c1) * (t2 - t1)
      t_total <- t_total + (t_cross - t1)
    } else if (!a1 && a2 && c2 != c1) {
      t_cross <- t1 + (mic - c1) / (c2 - c1) * (t2 - t1)
      t_total <- t_total + (t2 - t_cross)
    }
  }

  t_total
}

pk_model <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dA_depot <- -Ka * A_depot
    dA_central <- Ka * A_depot - (Cl / V1 + Q / V1) * A_central + (Q / V2) * A_periph
    dA_periph <- (Q / V1) * A_central - (Q / V2) * A_periph
    list(c(dA_depot, dA_central, dA_periph), Cp = A_central / V1)
  })
}

simulate_pk <- function(dose, interval, n_doses, pk_params = PK_PARAMS, dt = 0.25) {
  params <- c(
    Ka = pk_params$Ka_pop,
    Cl = pk_params$Cl_pop,
    V1 = pk_params$V1_pop,
    V2 = pk_params$V2_pop,
    Q = pk_params$Q_pop
  )

  t_end <- (n_doses - 1) * interval + 72
  times <- seq(0, t_end, by = dt)
  state <- c(A_depot = 0, A_central = 0, A_periph = 0)

  dose_times <- seq(0, (n_doses - 1) * interval, by = interval) + pk_params$Tlag_pop
  eventdat <- data.frame(
    var = rep("A_depot", n_doses),
    time = dose_times,
    value = rep(dose, n_doses),
    method = rep("add", n_doses)
  )

  out <- ode(
    y = state,
    times = times,
    func = pk_model,
    parms = params,
    events = list(data = eventdat),
    method = "lsoda"
  )

  result <- as.data.frame(out)
  result$Cp <- result$A_central / pk_params$V1_pop
  result$dose <- dose
  result
}

# =============================================================================
# UI
# =============================================================================

ui <- page_sidebar(
  title = "Modelo PK de florfenicol em tilápias (28 °C)",
  theme = bs_theme(
    bootswatch = "flatly",
    base_font = font_collection("Open Sans", "Helvetica Neue", "sans-serif"),
    heading_font = font_collection("Open Sans", "Helvetica Neue", "sans-serif")
  ),

  sidebar = sidebar(
    width = 340,

    h5("Protocolo Terapêutico"),
    sliderInput("dose", "Dose (mg/kg):", min = 5, max = 50, value = 15, step = 1),
    sliderInput("interval", "Intervalo entre doses (h):", min = 12, max = 48, value = 24, step = 6),
    sliderInput("n_doses", "Número de doses:", min = 1, max = 15, value = 10, step = 1),

    h5("Linha de Referência (CIM)"),
    numericInput("mic_line", "CIM (µg/mL, definida pelo usuário):",
                 value = 1, min = 0.001, max = 100, step = 0.1),
    checkboxInput("show_mic", "Mostrar linha CIM", value = TRUE),
    p("O artigo cita faixas de CIM, não um valor único: Streptococcus agalactiae 0.125-16 e Aeromonas spp. 0.125-4 µg/mL.",
      style = "font-size:0.78em; color:#666; margin-top:6px;"),

    hr(),

    tags$div(
      style = "background:#eafaf1; border-left:3px solid #27ae60; padding:8px 10px; border-radius:4px;",
      h6("Parâmetros do Modelo Publicado (28 °C)",
         style = "margin:0 0 4px 0; color:#1e8449;"),
      p("Valores publicados por Tameirão et al. (2022).",
        style = "font-size:0.78em; color:#555; margin:0 0 8px 0;"),
      numericInput("Tlag_pop", "Tlag (h):", value = 7.61, min = 0, max = 48, step = 0.1),
      numericInput("Ka_pop", "Ka (1/h):", value = 0.83, min = 0.001, max = 50, step = 0.01),
      numericInput("Cl_pop", "Cl (L/h/kg):", value = 0.0039, min = 0.0001, max = 10, step = 0.0001),
      numericInput("V1_pop", "V1 (L/kg):", value = 0.31, min = 0.001, max = 10, step = 0.01),
      numericInput("V2_pop", "V2 (L/kg):", value = 0.13, min = 0.001, max = 10, step = 0.01),
      numericInput("Q_pop", "Q (L/h/kg):", value = 0.057, min = 0, max = 10, step = 0.001)
    ),

    br(),
    actionButton("simulate", "Simular", class = "btn-primary btn-lg w-100")
  ),

  navset_card_tab(
    nav_panel(
      "Concentração plasmática",
      card_header("Concentração plasmática de florfenicol ao longo do tempo"),
      plotOutput("pk_plot", height = "500px"),
      card_footer(
        "Modelo publicado: oral, absorção de 1ª ordem com Tlag, 2 compartimentos e eliminação linear. ",
        "Temperatura fixa na referência de 28 °C."
      )
    ),

    nav_panel(
      "Comparação de doses",
      card_header("Comparação entre doses do artigo (10, 15, 20 mg/kg)"),
      plotOutput("dose_comparison_plot", height = "500px"),
      card_footer(
        "Simulação comparativa com 10, 15 e 20 mg/kg a cada 24 h por 10 dias, como descrito no artigo. ",
        "Paleta ajustada para maior contraste."
      )
    ),

    nav_panel(
      "Parâmetros PK",
      card_header("Parâmetros farmacocinéticos publicados"),
      layout_columns(
        col_widths = c(6, 6),
        tableOutput("pk_params_table"),
        tags$div(
          style = "padding: 8px 12px;",
          h5("Parâmetros"),
          tags$dl(
            tags$dt("Tlag"),
            tags$dd("Tempo de latência entre a administração e o início da absorção."),
            tags$dt("Ka"),
            tags$dd("Constante de absorção de primeira ordem."),
            tags$dt("Cl"),
            tags$dd("Clearance do compartimento central."),
            tags$dt("V1"),
            tags$dd("Volume aparente do compartimento central."),
            tags$dt("V2"),
            tags$dd("Volume aparente do compartimento periférico."),
            tags$dt("Q"),
            tags$dd("Clearance intercompartimental entre V1 e V2."),
            tags$dt("Temperatura de referência"),
            tags$dd("Temperatura usada como referência para os parâmetros publicados no artigo.")
          )
        )
      ),
      card_footer("Tabela com os parâmetros estruturais publicados para a referência de 28 °C.")
    ),

    nav_panel(
      "Métricas PK",
      card_header("Métricas farmacocinéticas"),
      tableOutput("pk_metrics_table"),
      card_footer(
        "Cmax e Cmin estimados no último intervalo. ",
        "AUC calculada por método trapezoidal. ",
        "T>CIM calculado por integração temporal com interpolação linear."
      )
    ),

    nav_panel(
      "Sobre",
      card_header("Sobre o modelo"),
      card_body(
        h4("Modelo farmacocinético de florfenicol em tilápias do Nilo"),
        p("Este aplicativo implementa o modelo estrutural descrito por:"),
        tags$blockquote(
          "Tameirão ER, Rubim FM, Felix LA, Gonzaga LWF, Brandão HM, ",
          "Murgas LDS, Ferrante M. (2022). Modelo farmacocinético de florfenicol ",
          "en tilapias (Oreochromis niloticus) sometidas a diferentes temperaturas ",
          "de crianza. Rev Inv Vet Peru, 33(6): e22433. ",
          tags$a("https://doi.org/10.15381/rivep.v33i6.22433",
                 href = "https://doi.org/10.15381/rivep.v33i6.22433",
                 target = "_blank")
        ),
        hr(),
        h5("CIM no artigo"),
        p("O artigo não define um único valor de CIM para usar no gráfico. Ele cita faixas de suscetibilidade por patógeno:"),
        tags$ul(
          tags$li(tags$em("Streptococcus agalactiae"), ": 0.125-16 µg/mL"),
          tags$li(tags$em("Aeromonas"), " spp.: 0.125-4 µg/mL")
        ),
        p("Por isso, o campo de CIM neste app é definido pelo usuário e usa 1 µg/mL apenas como valor inicial."),
        hr(),
        p("Nota: esta é uma ferramenta de apoio didático e de pesquisa. Decisões terapêuticas devem ser tomadas por médico veterinário qualificado.",
          style = "color:#c0392b; font-weight:bold;")
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  get_pk_params <- reactive({
    req(input$Tlag_pop, input$Ka_pop, input$Cl_pop,
        input$V1_pop, input$V2_pop, input$Q_pop)
    list(
      Tlag_pop = input$Tlag_pop,
      Ka_pop = input$Ka_pop,
      Cl_pop = input$Cl_pop,
      V1_pop = input$V1_pop,
      V2_pop = input$V2_pop,
      Q_pop = input$Q_pop,
      T_ref = 28
    )
  })

  sim_data <- eventReactive(input$simulate, {
    pk_params <- get_pk_params()
    validate_pk_params(pk_params)

    simulate_pk(
      dose = input$dose,
      interval = input$interval,
      n_doses = input$n_doses,
      pk_params = pk_params
    )
  }, ignoreNULL = FALSE)

  dose_comparison_data <- eventReactive(input$simulate, {
    pk_params <- get_pk_params()
    validate_pk_params(pk_params)

    doses <- c(10, 15, 20)

    bind_rows(lapply(doses, function(d) {
      res <- simulate_pk(
        dose = d,
        interval = 24,
        n_doses = 10,
        pk_params = pk_params
      )
      res$dose_label <- paste0(d, " mg/kg")
      res
    }))
  }, ignoreNULL = FALSE)

  metrics_data <- eventReactive(input$simulate, {
    pk_params <- get_pk_params()
    validate_pk_params(pk_params)

    mic <- if (!is.na(input$mic_line) && input$mic_line > 0) input$mic_line else 1
    res <- simulate_pk(
      dose = input$dose,
      interval = input$interval,
      n_doses = input$n_doses,
      pk_params = pk_params
    )

    t_last <- (input$n_doses - 1) * input$interval
    last_int <- res %>%
      filter(time >= t_last, time <= t_last + input$interval) %>%
      arrange(time)

    req(nrow(last_int) > 1)

    auc <- sum(diff(last_int$time) *
                 (head(last_int$Cp, -1) + tail(last_int$Cp, -1)) / 2)

    data.frame(
      `Dose (mg/kg)` = input$dose,
      `Intervalo (h)` = input$interval,
      `N doses` = input$n_doses,
      `Cmax (ug/mL)` = round(max(last_int$Cp, na.rm = TRUE), 2),
      `Cmin (ug/mL)` = round(min(last_int$Cp, na.rm = TRUE), 2),
      `AUC (ug.h/mL)` = round(auc, 2),
      `T>CIM (h)` = round(time_above_mic(last_int$time, last_int$Cp, mic), 1),
      check.names = FALSE
    )
  }, ignoreNULL = FALSE)

  output$pk_plot <- renderPlot({
    data <- sim_data()
    req(nrow(data) > 0)

    p <- ggplot(data, aes(x = time, y = Cp)) +
      geom_line(color = MAIN_CURVE_COLOR, linewidth = 1.2) +
      labs(
        x = "Tempo (h)",
        y = "Concentração Plasmática (µg/mL)",
        title = paste0("Florfenicol ", input$dose, " mg/kg, a cada ",
                       input$interval, " h, ", input$n_doses, " doses")
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    if (input$show_mic && !is.na(input$mic_line) && input$mic_line > 0) {
      label_y <- input$mic_line + 0.03 * max(data$Cp, na.rm = TRUE)

      p <- p +
        geom_hline(yintercept = input$mic_line,
                   linetype = "dashed", color = "#C62828", linewidth = 0.8) +
        annotate(
          "text",
          x = max(data$time) * 0.82,
          y = label_y,
          label = paste0("CIM = ", input$mic_line, " µg/mL"),
          color = "#C62828",
          size = 4
        )
    }

    p
  })

  output$dose_comparison_plot <- renderPlot({
    data <- dose_comparison_data()
    req(nrow(data) > 0)

    ggplot(data, aes(x = time, y = Cp, color = dose_label)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = DOSE_COLORS, name = "Dose") +
      labs(
        x = "Tempo (h)",
        y = "Concentração Plasmática (µg/mL)",
        title = "Comparação entre doses do artigo"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  })

  output$pk_params_table <- renderTable({
    pk_params <- get_pk_params()
    validate_pk_params(pk_params)

    data.frame(
      Parametro = c("Tlag", "Ka", "Cl", "V1", "V2", "Q", "Temperatura de referência"),
      Valor = c(
        round(pk_params$Tlag_pop, 2),
        round(pk_params$Ka_pop, 3),
        round(pk_params$Cl_pop, 5),
        round(pk_params$V1_pop, 2),
        round(pk_params$V2_pop, 2),
        round(pk_params$Q_pop, 3),
        pk_params$T_ref
      ),
      Unidade = c("h", "1/h", "L/h/kg", "L/kg", "L/kg", "L/h/kg", "C"),
      check.names = FALSE
    )
  }, digits = 4)

  output$pk_metrics_table <- renderTable({
    metrics_data()
  }, digits = 2)
}

# =============================================================================
shinyApp(ui, server)
