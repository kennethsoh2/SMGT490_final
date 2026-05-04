library(tidyverse)
library(gt)
library(dplyr)
library(data.table)
library(lubridate)
library(parallel)
library(lme4)
library(sabRmetrics)
library(xgboost)
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

coors_pitch_summary <- tibble(
  pitch_bucket = c("FASTBALL","CURVE","CUTTER","SWEEPER","SINKER","CHANGEUP","SLIDER","SPLIT"),
  
  velocity_change = c(-0.03, -0.03, 0.01, -0.07, -0.02, -0.06, 0.05, -0.10),
  
  spin_change = c(9.70, 20.40, 7.21, 7.27, 12.01, 17.63, 27.46, 5.02),
  
  ivb_change = c(-3.21, 1.95, -1.82, -1.60, -1.50, -1.20, -0.59, -0.26),
  
  hb_change = c(2.00, -1.75, 0.28, -2.61, 3.59, 3.59, -0.60, 3.40)
)

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

pitch_physics_raw <- pitch_metrics %>%
  mutate(
    velo = release_speed_mph,
    spin = release_spin_rate,
    ivb = induced_vert_break,
    hb = horz_break,
    outcome = delta_run_exp
  ) %>%
  filter(
    !is.na(outcome),
    !pitch_type %in% c("FA","CS","EP")
  ) %>%
  mutate(
    hb = ifelse(pitch_hand == "R", -hb, hb)
  ) %>%
  filter(
    velo > 70,
    velo < 105,
    abs(ivb) < 30,
    abs(hb) < 30
  )

pitch_physics_raw <- pitch_physics_raw %>%
  mutate(
    pitch_bucket = case_when(
      pitch_type == "FF" ~ "FASTBALL",
      pitch_type == "SI" ~ "SINKER",
      pitch_type == "FC" ~ "CUTTER",
      pitch_type == "SL" ~ "SLIDER",
      pitch_type %in% c("SV","ST") ~ "SWEEPER",
      pitch_type %in% c("CU","KC") ~ "CURVE",
      pitch_type == "CH" ~ "CHANGEUP",
      pitch_type %in% c("FS","FO") ~ "SPLIT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pitch_bucket))

pitch_physics_raw <- pitch_physics_raw %>%
  select(
    pitcher_id,
    pitch_hand,
    year,
    coors,
    pitch_bucket,
    velo,
    spin,
    ivb,
    hb,
    VAA,
    HAA,
    extension,
    arm_angle,
    outcome
  ) %>%
  drop_na()

saveRDS(pitch_physics_raw, "pitch_data.rds")

stuff_train <- pitch_physics_raw %>%
  filter(coors == 0) %>%
  mutate(pitch_bucket = as.factor(pitch_bucket))

pitch_levels <- levels(stuff_train$pitch_bucket)

feature_formula <- ~ velo + spin + ivb + hb + VAA + HAA + extension + arm_angle + pitch_bucket

saveRDS(pitch_levels, "pitch_levels.rds")
saveRDS(feature_formula, "feature_formula.rds")

X <- model.matrix(
  feature_formula,
  data = stuff_train
)[, -1]

y <- stuff_train$outcome

dtrain <- xgb.DMatrix(data = X, label = y)

params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.05,
  max_depth = 6,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 5,
  gamma = 0.1,
  lambda = 1,
  alpha = 0.5,
  nthread = parallel::detectCores()
)

stuff_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  verbose = 1
)

saveRDS(stuff_model, "stuff_model_final.rds")

predict_stuff_value <- function(data, model, pitch_levels) {
  
  data <- data %>%
    mutate(pitch_bucket = factor(pitch_bucket, levels = pitch_levels))
  
  X <- model.matrix(
    feature_formula,
    data = data
  )[, -1]
  
  dtest <- xgb.DMatrix(data = X)
  
  preds <- predict(model, dtest)
  
  return(preds)
}

train_preds <- predict_stuff_value(stuff_train, stuff_model, pitch_levels)

lower_bound <- quantile(train_preds, 0.01, na.rm = TRUE)
upper_bound <- quantile(train_preds, 0.99, na.rm = TRUE)

train_preds_clipped <- pmax(pmin(train_preds, upper_bound), lower_bound)

stuff_mean <- mean(train_preds_clipped, na.rm = TRUE)
stuff_sd   <- sd(train_preds_clipped, na.rm = TRUE)

saveRDS(list(
  mean = stuff_mean,
  sd = stuff_sd,
  lower = lower_bound,
  upper = upper_bound
), "stuff_scaling.rds")

compute_stuff_plus <- function(preds, mean_val, sd_val) {
  100 - 30 * (preds - mean_val) / sd_val
}

apply_coors_adjustment <- function(data, coors_table) {
  
  data %>%
    left_join(coors_table, by = "pitch_bucket") %>%
    mutate(
      velo = velo + velocity_change,
      spin = spin + spin_change,
      ivb  = ivb  + ivb_change,
      hb   = hb   + hb_change
    ) %>%
    select(-velocity_change, -spin_change, -ivb_change, -hb_change)
}

test_pitcher <- pitch_physics_raw %>%
  filter(pitcher_id == 657277, year == 2025)

baseline_preds <- predict_stuff_value(test_pitcher, stuff_model, pitch_levels)

baseline_preds <- pmax(pmin(baseline_preds, upper_bound), lower_bound)

baseline_stuff <- compute_stuff_plus(
  baseline_preds,
  stuff_mean,
  stuff_sd
)

coors_pitcher <- apply_coors_adjustment(test_pitcher, coors_pitch_summary)

coors_preds <- predict_stuff_value(coors_pitcher, stuff_model, pitch_levels)

coors_preds <- pmax(pmin(coors_preds, upper_bound), lower_bound)

coors_stuff <- compute_stuff_plus(
  coors_preds,
  stuff_mean,
  stuff_sd
)

results <- test_pitcher %>%
  mutate(
    baseline_stuff = baseline_stuff,
    coors_stuff = coors_stuff
  ) %>%
  group_by(pitch_bucket) %>%
  summarize(
    baseline = mean(baseline_stuff),
    coors = mean(coors_stuff),
    delta = coors - baseline,
    n = n()
  )

results

ggplot(results, aes(x = pitch_bucket)) +
  geom_col(aes(y = baseline, fill = "Baseline"), position = "dodge") +
  geom_col(aes(y = coors, fill = "Coors"), position = "dodge") +
  labs(
    title = "Stuff+ Change at Coors Field",
    y = "Stuff+",
    x = "Pitch Type"
  ) +
  theme_minimal()
