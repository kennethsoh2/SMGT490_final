library(shiny)
library(tidyverse)
library(xgboost)
library(ggplot2)

stuff_model <- readRDS("data/stuff_model_final.rds")
pitch_levels <- readRDS("data/pitch_levels.rds")
feature_formula <- readRDS("data/feature_formula.rds")

scaling <- readRDS("data/stuff_scaling.rds")
stuff_mean <- scaling$mean
stuff_sd <- scaling$sd
lower_bound <- scaling$lower
upper_bound <- scaling$upper

pitch_physics_raw <- readRDS("data/pitch_data_app.rds")

pitch_physics_raw <- pitch_physics_raw %>%
  mutate(pitcher_id = as.character(pitcher_id))

player_lookup_raw <- read_csv("data/player.csv")

pitcher_lookup <- player_lookup_raw %>%
  transmute(
    pitcher_id = as.character(player_id),
    pitcher_name = name_full
  ) %>%
  distinct() %>%
  drop_na() %>%
  arrange(pitcher_name)

pitcher_lookup <- pitcher_lookup %>%
  filter(pitcher_id %in% pitch_physics_raw$pitcher_id)

coors_pitch_summary <- tibble(
  pitch_bucket = c("FASTBALL","CURVE","CUTTER","SWEEPER","SINKER","CHANGEUP","SLIDER","SPLIT"),
  velocity_change = c(-0.03, -0.03, 0.01, -0.07, -0.02, -0.06, 0.05, -0.10),
  spin_change = c(9.70, 20.40, 7.21, 7.27, 12.01, 17.63, 27.46, 5.02),
  ivb_change = c(-3.21, 1.95, -1.82, -1.60, -1.50, -1.20, -0.59, -0.26),
  hb_change = c(2.00, -1.75, 0.28, -2.61, 3.59, 3.59, -0.60, 3.40)
)

predict_stuff_value <- function(data, model, pitch_levels) {
  data <- data %>%
    mutate(pitch_bucket = factor(pitch_bucket, levels = pitch_levels))
  
  X <- model.matrix(feature_formula, data = data)[, -1]
  dtest <- xgb.DMatrix(data = X)
  predict(model, dtest)
}

compute_stuff_plus <- function(preds) {
  preds <- pmax(pmin(preds, upper_bound), lower_bound)
  100 - 30 * (preds - stuff_mean) / stuff_sd
}

apply_coors_adjustment <- function(data) {
  data %>%
    left_join(coors_pitch_summary, by = "pitch_bucket") %>%
    mutate(
      velo = velo + velocity_change,
      spin = spin + spin_change,
      ivb  = ivb  + ivb_change,
      hb   = hb   + hb_change
    ) %>%
    select(-velocity_change, -spin_change, -ivb_change, -hb_change)
}

ui <- fluidPage(
  
  titlePanel("Coors Field Pitch Profile Simulator"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(
        "pitcher",
        "Select Pitcher",
        choices = setNames(
          pitcher_lookup$pitcher_id,
          pitcher_lookup$pitcher_name
        )
      ),
      
      selectInput(
        "year",
        "Select Year",
        choices = sort(unique(pitch_physics_raw$year))
      )
    ),
    
    mainPanel(
      
      h3(textOutput("player_name")),
      
      h4("Baseline Pitch Metrics (Season Averages)"),
      tableOutput("baseline_table"),
      
      h4("Coors-Adjusted Pitch Metrics (Simulated)"),
      tableOutput("coors_table"),
      
      h4("Pitch Movement (Baseline)"),
      plotOutput("baseline_plot", height = "400px"),
      
      h4("Pitch Movement (Coors Adjusted)"),
      plotOutput("coors_plot", height = "400px")
    )
  )
)

server <- function(input, output) {
  
  output$player_name <- renderText({
    name <- pitcher_lookup %>%
      filter(pitcher_id == input$pitcher) %>%
      pull(pitcher_name) %>% 
      first()
    
    paste("Pitcher:", name, "| Year:", input$year)
  })
  
  pitcher_data <- reactive({
    pitch_physics_raw %>%
      filter(
        pitcher_id == input$pitcher,
        year == input$year
      )
  })
  
  analysis_data <- reactive({
    
    df <- pitcher_data()
    if (nrow(df) == 0) return(NULL)
  
    base_preds <- predict_stuff_value(df, stuff_model, pitch_levels)
    base_stuff <- compute_stuff_plus(base_preds)
    
    baseline_df <- df %>%
      mutate(stuff = base_stuff)
    
    coors_df <- apply_coors_adjustment(df)
    coors_preds <- predict_stuff_value(coors_df, stuff_model, pitch_levels)
    coors_stuff <- compute_stuff_plus(coors_preds)
    
    coors_df <- coors_df %>%
      mutate(stuff = coors_stuff)
    
    list(
      baseline = baseline_df,
      coors = coors_df
    )
  })
  
  clean_pitch_names <- function(x) {
    recode(x,
           "FASTBALL" = "Fastball",
           "SINKER" = "Sinker",
           "CUTTER" = "Cutter",
           "SLIDER" = "Slider",
           "SWEEPER" = "Sweeper",
           "CURVE" = "Curveball",
           "CHANGEUP" = "Changeup",
           "SPLIT" = "Splitter"
    )
  }
  
  format_1_decimal <- function(x) {
    format(round(x, 1), nsmall = 1)
  }
  
  format_int <- function(x) {
    format(round(x), trim = TRUE)
  }
  
  output$baseline_table <- renderTable({
    
    data <- analysis_data()
    if (is.null(data)) return(NULL)
    
    df <- data$baseline
    total_pitches <- nrow(df)
    
    df %>%
      group_by(pitch_bucket) %>%
      summarize(
        Usage = paste0(round(100 * n() / total_pitches), "%"),
        Velo = format_1_decimal(mean(velo, na.rm = TRUE)),
        Spin = format_int(mean(spin, na.rm = TRUE)),
        IVB = format_1_decimal(mean(ivb, na.rm = TRUE)),
        HB = format_1_decimal(mean(hb, na.rm = TRUE)),
        `Stuff+` = format_int(mean(stuff, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        Pitch = clean_pitch_names(pitch_bucket)
      ) %>%
      select(Pitch, Usage, Velo, Spin, IVB, HB, `Stuff+`) %>%
      arrange(desc(as.numeric(gsub("%", "", Usage))))
  })
  
  output$coors_table <- renderTable({
    
    data <- analysis_data()
    if (is.null(data)) return(NULL)
    
    df <- data$coors
    total_pitches <- nrow(df)
    
    df %>%
      group_by(pitch_bucket) %>%
      summarize(
        Usage = paste0(round(100 * n() / total_pitches), "%"),
        Velo = format_1_decimal(mean(velo, na.rm = TRUE)),
        Spin = format_int(mean(spin, na.rm = TRUE)),
        IVB = format_1_decimal(mean(ivb, na.rm = TRUE)),
        HB = format_1_decimal(mean(hb, na.rm = TRUE)),
        `Stuff+` = format_int(mean(stuff, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        Pitch = clean_pitch_names(pitch_bucket)
      ) %>%
      select(Pitch, Usage, Velo, Spin, IVB, HB, `Stuff+`) %>%
      arrange(desc(as.numeric(gsub("%", "", Usage))))
  })
  
  output$baseline_plot <- renderPlot({
    
    data <- analysis_data()
    if (is.null(data)) return(NULL)
    
    plot_df <- data$baseline %>%
      mutate(
        hb_plot = ifelse(!is.na(pitch_hand) & pitch_hand == "L", -hb, hb)
      )
    
    if (nrow(plot_df) == 0) return(NULL)
    
    ggplot(plot_df, aes(x = hb_plot, y = ivb, color = pitch_bucket)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      coord_fixed() +
      labs(
        title = "Pitch Movement (Baseline)",
        x = "Horizontal Break",
        y = "Induced Vertical Break"
      ) +
      theme_minimal()
  })
  
  output$coors_plot <- renderPlot({
    
    data <- analysis_data()
    if (is.null(data)) return(NULL)
    
    plot_df <- data$coors %>%
      mutate(
        hb_plot = ifelse(!is.na(pitch_hand) & pitch_hand == "L", -hb, hb)
      )
    
    if (nrow(plot_df) == 0) return(NULL)
    
    ggplot(plot_df, aes(x = hb_plot, y = ivb, color = pitch_bucket)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      coord_fixed() +
      labs(
        title = "Pitch Movement (Coors Adjusted)",
        x = "Horizontal Break",
        y = "Induced Vertical Break"
      ) +
      theme_minimal()
  })
}

shinyApp(ui, server)
