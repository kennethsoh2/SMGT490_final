library(tidyverse)
library(gt)
library(dplyr)
library(data.table)
library(lubridate)
library(parallel)
library(lme4)
library(mclust)
library(sabRmetrics)
library(ggplot2)
library(broom.mixed)
library(scales)
library(here)

data_path <- here("data", "savant_data")

files <- list.files(
  path = data_path,
  pattern = "data_.*\\.csv",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No savant data files found in data/savant_data")
}

statcast <- rbindlist(
  lapply(files, function(f) {
    fread(f, showProgress = TRUE)
  }),
  use.names = TRUE,
  fill = TRUE
)

statcast <- statcast %>%
  mutate(
    game_date = as.Date(game_date),
    year = as.integer(year),
    coors = ifelse(home_team == "COL", 1, 0),
    rockies = case_when(
      home_team == "COL" & inning_topbot == "Bot" ~ 1,
      away_team == "COL" & inning_topbot == "Top" ~ 1,
      TRUE ~ 0
    )
  ) %>%
  filter(game_type == "R")

pitch_metrics <- statcast %>%
  sabRmetrics::get_quadratic_coef(source = "baseballsavant") %>%
  sabRmetrics::get_trackman_metrics()

pitch_metrics <- pitch_metrics %>%
  mutate(
    release_speed_mph = release_speed / 1.4667,
    
    vx_plate = vx0 + ax * plate_time,
    vy_plate = vy0 + ay * plate_time,
    vz_plate = vz0 + az * plate_time,
    
    VAA = -atan(vz_plate / vy_plate) * 180 / pi,
    HAA = -atan(vx_plate / vy_plate) * 180 / pi
  ) %>%
  select(
    game_id,
    year,
    coors,
    rockies,
    pitcher_id,
    batter_id,
    pitch_hand,
    pitch_type,
    release_speed_mph,
    release_spin_rate,
    induced_vert_break,
    horz_break,
    extension,
    arm_angle,
    VAA,
    HAA,
    delta_run_exp
    
  )

pitch_physics <- pitch_metrics %>%
  mutate(
    velo = release_speed_mph,
    spin = release_spin_rate,
    ivb = induced_vert_break,
    hb = horz_break,
    extension = extension,
    arm_angle = arm_angle,
    outcome = delta_run_exp
  ) %>%
  filter(
    !is.na(outcome),
    !pitch_type %in% c("FA","CS","EP")
  )

pitch_physics <- pitch_physics %>%
  mutate(
    hb = ifelse(pitch_hand == "R", -hb, hb)
  )

pitch_physics <- pitch_physics %>%
  filter(
    velo > 70,
    velo < 105,
    abs(ivb) < 30,
    abs(hb) < 30
  )

pitch_physics <- pitch_physics %>%
  mutate(
    velo_raw = velo,
    spin_raw = spin,
    ivb_raw  = ivb,
    hb_raw   = hb
  ) %>%
  mutate(
    velo = as.numeric(scale(velo_raw)),
    spin = as.numeric(scale(spin_raw)),
    ivb  = as.numeric(scale(ivb_raw)),
    hb   = as.numeric(scale(hb_raw))
  )

pitch_physics <- pitch_physics %>%
  mutate(
    pitch_bucket = case_when(
      pitch_type == "FF" ~ "FF",
      pitch_type == "SI" ~ "SI",
      pitch_type == "FC" ~ "FC",
      pitch_type == "SL" ~ "SL",
      pitch_type %in% c("SV","ST") ~ "SWEEPER",
      pitch_type %in% c("CU","KC") ~ "CURVE",
      pitch_type == "CH" ~ "CH",
      pitch_type %in% c("FS","FO") ~ "SPLIT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pitch_bucket))

set.seed(42)

pitch_sample <- pitch_physics %>%
  group_by(pitch_bucket) %>%
  sample_n(min(35000, n())) %>%
  ungroup()

pitch_matrix_sample <- pitch_sample %>%
  select(velo, spin, ivb, hb)

pitch_gmm <- Mclust(
  pitch_matrix_sample,
  G = 5:7
)

pitch_matrix_full <- pitch_physics %>%
  select(velo, spin, ivb, hb) %>%
  mutate(across(everything(),
                ~ifelse(is.infinite(.), NA, .)))

valid_idx <- complete.cases(pitch_matrix_full)

pitch_physics$shape_cluster <- NA

pitch_physics$shape_cluster[valid_idx] <- predict(
  pitch_gmm,
  newdata = pitch_matrix_full[valid_idx, ]
)$classification

pitch_physics$shape_cluster <- factor(pitch_physics$shape_cluster)

pitch_cluster_centers <- pitch_physics %>%
  filter(!is.na(shape_cluster)) %>%
  group_by(shape_cluster) %>%
  summarise(
    velo = mean(velo_raw),
    spin = mean(spin_raw),
    ivb = mean(ivb_raw),
    hb = mean(hb_raw),
    n = n(),
    .groups = "drop"
  )

print(pitch_cluster_centers)

pitch_cluster_table <- pitch_cluster_centers %>%
  mutate(
    velo = round(velo, 1),
    spin = round(spin, 0),
    ivb = round(ivb, 1),
    hb = round(hb, 1),
    n = format(n, big.mark = ",")
  ) %>%
  gt() %>%
  
  cols_label(
    shape_cluster = "Cluster",
    velo = "Velocity (mph)",
    spin = "Spin Rate (rpm)",
    ivb = "IVB (in)",
    hb = "HB (in)",
    n = "Count"
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  tab_header(
    title = "Pitch Shape Clusters",
    subtitle = "Cluster Centers Based on Velocity, Spin, and Movement"
  ) %>%
  
  fmt_number(
    columns = c(velo, spin, ivb, hb),
    decimals = 1
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  
  tab_options(
    table.font.size = 14,
    data_row.padding = px(6)
  )

pitch_cluster_table

gtsave(pitch_cluster_table, "pitch_shape_clusters.png")

ggplot(
  pitch_physics %>%
    filter(!is.na(shape_cluster)) %>%
    sample_n(40000),
  aes(x = hb, y = ivb, color = shape_cluster)
) +
  geom_point(alpha = 0.4, size = 1.2) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Pitch Shape Clusters (HB vs IVB)",
    x = "Horizontal Break",
    y = "Induced Vertical Break",
    color = "Cluster"
  )

pitch_mix <- pitch_physics %>%
  filter(!is.na(shape_cluster)) %>%
  group_by(pitcher_id, year, shape_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pitcher_id, year) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

pitch_mix_wide <- pitch_mix %>%
  select(pitcher_id, year, shape_cluster, pct) %>%
  pivot_wider(
    names_from = shape_cluster,
    values_from = pct,
    values_fill = 0,
    names_prefix = "cluster_"
  )

pitch_entropy <- pitch_mix %>%
  group_by(pitcher_id, year) %>%
  summarise(
    entropy = -sum(pct * log(pct + 1e-8)),
    .groups = "drop"
  )

pitch_features <- pitch_mix_wide %>%
  left_join(pitch_entropy, by = c("pitcher_id","year"))

pitch_counts <- pitch_physics %>%
  group_by(pitcher_id, year) %>%
  summarise(n_pitches = n(), .groups = "drop")

pitch_features <- pitch_features %>%
  left_join(pitch_counts, by = c("pitcher_id","year")) %>%
  filter(n_pitches >= 500)

pitch_matrix <- pitch_features %>%
  select(starts_with("cluster_"), entropy) %>%
  mutate(across(where(is.numeric),
                ~ifelse(is.infinite(.), NA, .))) %>%
  drop_na()

gmm_pitchers <- Mclust(
  pitch_matrix,
  G = 3:8
)

pitcher_clusters <- pitch_features %>%
  select(pitcher_id, year) %>%
  distinct() %>%
  mutate(cluster = factor(gmm_pitchers$classification))

pitcher_cluster_summary <- pitch_features %>%
  left_join(pitcher_clusters, by = c("pitcher_id","year")) %>%
  group_by(cluster) %>%
  summarise(
    across(starts_with("cluster_"), mean),
    entropy = mean(entropy),
    n_pitchers = n(),
    .groups = "drop"
  )

print(pitcher_cluster_summary)

pa_data <- statcast %>%
  group_by(game_id, at_bat_number) %>%
  summarise(
    runs = sum(delta_run_exp, na.rm = TRUE),
    pitcher_id = first(pitcher_id),
    batter_id = first(batter_id),
    year = first(year),
    coors = first(coors),
    rockies = first(rockies),
    .groups = "drop"
  )

pitcher_clusters <- pitcher_clusters %>%
  mutate(pitcher_id = as.character(pitcher_id))

pa_data_pitcher <- pa_data %>%
  mutate(pitcher_id = as.character(pitcher_id)) %>%
  inner_join(pitcher_clusters, by = c("pitcher_id","year"))

pitcher_cluster_model <- lmer(
  runs ~ coors * cluster + rockies + factor(year) +
    (1 | pitcher_id) +
    (1 | batter_id),
  data = pa_data_pitcher,
  REML = FALSE
)

summary(pitcher_cluster_model)

coef_df <- broom.mixed::tidy(pitcher_cluster_model, effects = "fixed") %>%
  
  filter(term != "(Intercept)") %>%
  
  mutate(
    term_clean = case_when(
      term == "coors" ~ "Coors × Cluster 1",
      
      term == "coors:cluster2" ~ "Coors × Cluster 2",
      term == "coors:cluster3" ~ "Coors × Cluster 3",
      term == "coors:cluster4" ~ "Coors × Cluster 4",
      term == "coors:cluster5" ~ "Coors × Cluster 5",
      term == "coors:cluster6" ~ "Coors × Cluster 6",
      term == "coors:cluster7" ~ "Coors × Cluster 7",
      term == "coors:cluster8" ~ "Coors × Cluster 8",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(!is.na(term_clean)) %>%
  
  mutate(
    term_clean = factor(term_clean, levels = c(
      "Coors × Cluster 1",
      "Coors × Cluster 2",
      "Coors × Cluster 3",
      "Coors × Cluster 4",
      "Coors × Cluster 5",
      "Coors × Cluster 6",
      "Coors × Cluster 7",
      "Coors × Cluster 8"
    )),
    
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2)
  ) %>%
  
  arrange(term_clean) %>%
  
  select(term_clean, estimate, std.error, statistic)

pitcher_cluster_table <- coef_df %>%
  gt() %>%
  
  cols_label(
    term_clean = "Variable",
    estimate = "Coefficient",
    std.error = "Std. Error",
    statistic = "t-value"
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  tab_header(
    title = "Pitcher Archetype Effects at Coors Field",
    subtitle = "Mixed Effects Model: Run Value per PA"
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  
  tab_options(
    table.font.size = 14,
    data_row.padding = px(6)
  )

pitcher_cluster_table

gtsave(pitcher_cluster_table, "pitcher_cluster_effects.png")

pitch_usage <- pitch_mix %>%
  select(pitcher_id, year, shape_cluster, pct) %>%
  pivot_wider(
    names_from = shape_cluster,
    values_from = pct,
    values_fill = 0,
    names_prefix = "pitch_"
  )

pitch_usage <- pitch_usage %>%
  mutate(pitcher_id = as.character(pitcher_id))

pa_data_usage <- pa_data %>%
  mutate(pitcher_id = as.character(pitcher_id)) %>%
  inner_join(pitch_usage, by = c("pitcher_id","year"))

pa_data_usage <- pa_data_usage %>%
  select(-pitch_1)

pitch_trait_model <- lmer(
  runs ~ coors * (
    pitch_2 + pitch_3 + pitch_4 +
      pitch_5 + pitch_6 + pitch_7
  ) + 
    rockies +
    factor(year) +
    (1 | pitcher_id) +
    (1 | batter_id),
  data = pa_data_usage,
  REML = FALSE
)

summary(pitch_trait_model)

coef_df <- broom.mixed::tidy(pitch_trait_model, effects = "fixed") %>%
  
  filter(term != "(Intercept)") %>%
  
  mutate(
    term_clean = case_when(
      term == "coors" ~ "Coors × Cluster 1",
      
      term == "coors:pitch_2" ~ "Coors × Cluster 2",
      term == "coors:pitch_3" ~ "Coors × Cluster 3",
      term == "coors:pitch_4" ~ "Coors × Cluster 4",
      term == "coors:pitch_5" ~ "Coors × Cluster 5",
      term == "coors:pitch_6" ~ "Coors × Cluster 6",
      term == "coors:pitch_7" ~ "Coors × Cluster 7",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(!is.na(term_clean)) %>%
  
  mutate(
    term_clean = factor(term_clean, levels = c(
      "Coors × Cluster 1",
      "Coors × Cluster 2",
      "Coors × Cluster 3",
      "Coors × Cluster 4",
      "Coors × Cluster 5",
      "Coors × Cluster 6",
      "Coors × Cluster 7"
    )),
    
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2)
  ) %>%
  
  arrange(term_clean) %>%
  
  select(term_clean, estimate, std.error, statistic)

pitch_trait_table <- coef_df %>%
  gt() %>%
  
  cols_label(
    term_clean = "Variable",
    estimate = "Coefficient",
    std.error = "Std. Error",
    statistic = "t-value"
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  tab_header(
    title = "Pitch Type Effects at Coors Field",
    subtitle = "Mixed Effects Model: Run Value per PA"
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  
  tab_options(
    table.font.size = 14,
    data_row.padding = px(6)
  )

pitch_trait_table

gtsave(pitch_trait_table, "pitch_trait_effects.png")

pitch_cols <- c(
  "cluster_1","cluster_2","cluster_3","cluster_4",
  "cluster_5","cluster_6","cluster_7"
)

midpoint <- 1/7

pitcher_cluster_summary_ordered <- pitcher_cluster_summary %>%
  select(
    cluster,
    cluster_2,
    cluster_7,
    cluster_1,
    cluster_6,
    cluster_4,
    cluster_3,
    cluster_5,
    entropy,
    n_pitchers
  )

pitcher_cluster_table <- pitcher_cluster_summary_ordered %>%
  gt() %>%

  cols_label(
    cluster = "Cluster",
    cluster_1 = "Slider/Cutter",
    cluster_2 = "4-Seam",
    cluster_3 = "Changeup",
    cluster_4 = "Sweeper",
    cluster_5 = "Splitter",
    cluster_6 = "Curveball",
    cluster_7 = "Sinker",
    entropy = "Entropy",
    n_pitchers = "# Pitchers"
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  fmt_percent(
    columns = all_of(pitch_cols),
    decimals = 0
  ) %>%

  data_color(
    columns = all_of(pitch_cols),
    colors = function(x) {
      scales::col_numeric(
        palette = c("#2C7BB6", "white", "#D7191C"),  # nicer blue/red
        domain = c(
          -max(abs(x - midpoint), na.rm = TRUE),
          max(abs(x - midpoint), na.rm = TRUE)
        )
      )(x - midpoint)
    }
  ) %>%

  fmt_number(
    columns = entropy,
    decimals = 2
  ) %>%

  fmt_number(
    columns = n_pitchers,
    decimals = 0,
    use_seps = TRUE
  ) %>%

  tab_header(
    title = "Pitcher Archetype Clusters",
    subtitle = "Average Pitch Mix (%) and Arsenal Diversity"
  ) %>%

  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%

  tab_options(
    table.font.size = 14,
    data_row.padding = px(6)
  )

pitcher_cluster_table

gtsave(pitcher_cluster_table, "pitcher_clusters_heatmap.png")
