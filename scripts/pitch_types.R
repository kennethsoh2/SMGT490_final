library(tidyverse)
library(gt)
library(dplyr)
library(data.table)
library(lubridate)
library(parallel)
library(scales)
library(lme4)
library(sabRmetrics)
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
    velo = as.numeric(scale(velo)),
    spin = as.numeric(scale(spin)),
    ivb = as.numeric(scale(ivb)),
    hb = as.numeric(scale(hb)),
    VAA = as.numeric(scale(VAA)),
    HAA = as.numeric(scale(HAA)),
    extension = as.numeric(scale(extension)),
    arm_angle = as.numeric(scale(arm_angle))
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

run_pitch_model <- function(data) {
  lmer(
    outcome ~ coors *
      (velo +
         spin +
          ivb +
          hb +
          VAA +
          HAA +
          extension +
          arm_angle
      ) +
      factor(year) +
      (1 | pitcher_id) +
      (1 | batter_id),
    data = data,
    REML = FALSE
  )
}

four_seam_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "FF")
)

sinker_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "SI")
)

cutter_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "FC")
)

slider_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "SL")
)

sweeper_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "SWEEPER")
)

curve_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "CURVE")
)

changeup_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "CH")
)

split_model <- run_pitch_model(
  pitch_physics %>% filter(pitch_bucket == "SPLIT")
)

summary(four_seam_model)
summary(sinker_model)
summary(cutter_model)
summary(slider_model)
summary(sweeper_model)
summary(curve_model)
summary(changeup_model)
summary(split_model)

extract_coors_effects <- function(model, pitch_name) {
  
  coef_table <- summary(model)$coefficients
  
  vars <- c("velo","spin","ivb","hb","VAA","HAA","extension","arm_angle")
  
  results <- lapply(vars, function(v) {
    
    base_est <- if (v %in% rownames(coef_table)) {
      coef_table[v, "Estimate"]
    } else NA
    
    inter_term <- paste0("coors:", v)
    
    inter_est <- if (inter_term %in% rownames(coef_table)) {
      coef_table[inter_term, "Estimate"]
    } else NA
    
    t_val <- if (inter_term %in% rownames(coef_table)) {
      coef_table[inter_term, "t value"]
    } else NA
    
    coors_effect <- base_est + inter_est
    
    tibble(
      variable = v,
      coors_effect = coors_effect,
      t_value = t_val
    )
  }) %>% bind_rows()
  
  wide <- results %>%
    pivot_wider(
      names_from = variable,
      values_from = c(coors_effect, t_value),
      names_glue = "{variable}_{.value}"
    )
  
  wide$pitch_type <- pitch_name
  
  wide
}

results_df <- bind_rows(
  extract_coors_effects(four_seam_model, "FF"),
  extract_coors_effects(sinker_model, "SI"),
  extract_coors_effects(cutter_model, "FC"),
  extract_coors_effects(slider_model, "SL"),
  extract_coors_effects(sweeper_model, "SWEEPER"),
  extract_coors_effects(curve_model, "CURVE"),
  extract_coors_effects(changeup_model, "CH"),
  extract_coors_effects(split_model, "SPLIT")
)

results_df <- results_df %>%
  select(
    pitch_type,
    
    velo_coors_effect, velo_t_value,
    spin_coors_effect, spin_t_value,
    ivb_coors_effect, ivb_t_value,
    hb_coors_effect, hb_t_value,
    VAA_coors_effect, VAA_t_value,
    HAA_coors_effect, HAA_t_value,
    extension_coors_effect, extension_t_value,
    arm_angle_coors_effect, arm_angle_t_value
  )

write.csv(results_df, "pitch_type_results.csv")

season_usage <- tibble(
  pitch_type = c("FF","SI","FC","SL","SWEEPER","CURVE","CH","SPLIT"),
  pitches_season = c(1200,900,700,800,600,500,600,400)
)

season_results <- results_df %>%
  left_join(season_usage, by = "pitch_type") %>%
  mutate(
    across(
      ends_with("coors_effect"),
      ~ . * pitches_season,
      .names = "{.col}_season_runs"
    )
  )

season_results <- season_results %>%
  mutate(across(-pitch_type, ~round(., 1)))

season_results <- season_results %>% 
  select(
    pitch_type,
    velo_coors_effect_season_runs,
    spin_coors_effect_season_runs,
    ivb_coors_effect_season_runs,
    hb_coors_effect_season_runs)

season_results <- season_results %>%
  mutate(
    pitch_type = recode(
      pitch_type,
      "FF" = "Four-Seam",
      "SI" = "Sinker",
      "FC" = "Cutter",
      "SL" = "Slider",
      "SWEEPER" = "Sweeper",
      "CURVE" = "Curveball",
      "CH" = "Changeup",
      "SPLIT" = "Splitter"
    )
  )

season_results <- season_results %>%
  gt() %>%
  
  cols_label(
    pitch_type = "Pitch Type",
    velo_coors_effect_season_runs = "Velocity",
    spin_coors_effect_season_runs = "Spin",
    ivb_coors_effect_season_runs = "IVB",
    hb_coors_effect_season_runs = "HB"
  ) %>%
  
  tab_header(
    title = "Seasonal Run Impact of Pitch Traits at Coors Field",
    subtitle = "Runs Saved/Lost Over a Season From a 1 SD Change in Each Pitch Characteristic"
  ) %>%
  
  cols_align(
    align = "center",
    -pitch_type
  ) %>%
  
  data_color(
    columns = -pitch_type,
    colors = scales::col_numeric(
      palette = c("#b2182b", "white", "#2166ac"),
      domain = NULL
    )
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

season_results
gtsave(season_results, "season_results.png")

pitch_eda <- pitch_metrics %>%
  mutate(
    velo = release_speed_mph,
    spin = release_spin_rate,
    ivb = induced_vert_break,
    hb = horz_break
  ) %>%
  filter(
    velo > 70,
    velo < 105,
    abs(ivb) < 30,
    abs(hb) < 30
  ) %>%
  mutate(
    hb = ifelse(pitch_hand == "L", -hb, hb),
    
    pitch_bucket = case_when(
      pitch_type == "FF" ~ "Four-Seam",
      pitch_type == "SI" ~ "Sinker",
      pitch_type == "FC" ~ "Cutter",
      pitch_type == "SL" ~ "Slider",
      pitch_type %in% c("SV","ST") ~ "Sweeper",
      pitch_type %in% c("CU","KC") ~ "Curveball",
      pitch_type == "CH" ~ "Changeup",
      pitch_type %in% c("FS","FO") ~ "Splitter",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pitch_bucket))


pitcher_env_compare <- pitch_eda %>%
  group_by(pitch_bucket, pitcher_id, year) %>%
  summarise(
    
    n_noncoors = sum(coors == 0),
    n_coors = sum(coors == 1),
    
    velo_noncoors = mean(velo[coors == 0], na.rm = TRUE),
    velo_coors = mean(velo[coors == 1], na.rm = TRUE),
    
    spin_noncoors = mean(spin[coors == 0], na.rm = TRUE),
    spin_coors = mean(spin[coors == 1], na.rm = TRUE),
    
    ivb_noncoors = mean(ivb[coors == 0], na.rm = TRUE),
    ivb_coors = mean(ivb[coors == 1], na.rm = TRUE),
    
    hb_noncoors = mean(hb[coors == 0], na.rm = TRUE),
    hb_coors = mean(hb[coors == 1], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  filter(
    n_noncoors >= 20,
    n_coors >= 5
  ) %>%
  mutate(
    velo_diff = velo_coors - velo_noncoors,
    spin_diff = spin_coors - spin_noncoors,
    ivb_diff = ivb_coors - ivb_noncoors,
    hb_diff = hb_coors - hb_noncoors
  )

coors_pitch_summary <- pitcher_env_compare %>%
  group_by(pitch_bucket) %>%
  summarise(
    
    velocity_change = mean(velo_diff, na.rm = TRUE),
    spin_change = mean(spin_diff, na.rm = TRUE),
    ivb_change = mean(ivb_diff, na.rm = TRUE),
    hb_change = mean(hb_diff, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(abs(ivb_change)))


coors_pitch_summary <- coors_pitch_summary %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  gt() %>%
  cols_label(
    pitch_bucket = "Pitch Type",
    velocity_change = "Δ Velocity (mph)",
    spin_change = "Δ Spin (rpm)",
    ivb_change = "Δ IVB (inches)",
    hb_change = "Δ HB (inches)"
  ) %>%
  tab_header(
    title = "Pitch Changes at Coors Field",
    subtitle = "Average Difference Between Coors and Non-Coors Environments"
  )

coors_pitch_summary
gtsave(coors_pitch_summary, "pitch_changes.png")