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

setwd("~/R/SMGT 490")

# --------------------------------------------------
# LOAD STATCAST DATA
# --------------------------------------------------

data_path <- "savant_data"

files <- list.files(
  path = data_path,
  pattern = "data_.*\\.csv",
  full.names = TRUE
)

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

# --------------------------------------------------
# RECONSTRUCT PITCH TRAJECTORY
# --------------------------------------------------

pitch_metrics <- statcast %>%
  sabRmetrics::get_quadratic_coef(source = "baseballsavant") %>%
  sabRmetrics::get_trackman_metrics()

# --------------------------------------------------
# COMPUTE APPROACH ANGLES
# --------------------------------------------------

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
    
    # identifiers
    game_id,
    year,
    
    # park indicator
    coors,
    rockies,
    
    # pitcher / batter info
    pitcher_id,
    batter_id,
    pitch_hand,
    pitch_type,
    
    # pitch physics
    release_speed_mph,
    release_spin_rate,
    induced_vert_break,
    horz_break,
    extension,
    arm_angle,
    
    # approach angles
    VAA,
    HAA,
    
    # modeling outcome
    delta_run_exp
    
  )

# --------------------------------------------------
# BUILD PITCH PHYSICS DATASET
# --------------------------------------------------

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

# --------------------------------------------------
# ACCOUNT FOR HANDEDNESS (HB SIGN)
# --------------------------------------------------

pitch_physics <- pitch_physics %>%
  mutate(
    hb = ifelse(pitch_hand == "R", -hb, hb)
  )

# --------------------------------------------------
# REMOVE EXTREME OUTLIERS
# --------------------------------------------------

pitch_physics <- pitch_physics %>%
  filter(
    velo > 70,
    velo < 105,
    abs(ivb) < 30,
    abs(hb) < 30,
  )

# --------------------------------------------------
# STANDARDIZE PHYSICS VARIABLES
# --------------------------------------------------

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

# --------------------------------------------------
# CREATE PITCH BUCKETS
# --------------------------------------------------

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

# --------------------------------------------------
# PITCH SHAPE CLUSTERING (STRATIFIED SAMPLE)
# --------------------------------------------------

set.seed(42)

# Stratified sampling by pitch_bucket
pitch_sample <- pitch_physics %>%
  group_by(pitch_bucket) %>%
  sample_n(min(35000, n())) %>%
  ungroup()

# Select clustering features
pitch_matrix_sample <- pitch_sample %>%
  select(velo, spin, ivb, hb)

# Run GMM on sampled data
pitch_gmm <- Mclust(
  pitch_matrix_sample,
  G = 5:8
)

# --------------------------------------------------
# ASSIGN CLUSTERS TO FULL DATASET (NA-SAFE)
# --------------------------------------------------

pitch_matrix_full <- pitch_physics %>%
  select(velo, spin, ivb, hb) %>%
  mutate(across(everything(),
                ~ifelse(is.infinite(.), NA, .)))

# Identify valid rows (no NA)
valid_idx <- complete.cases(pitch_matrix_full)

# Initialize cluster column
pitch_physics$shape_cluster <- NA

# Predict only on valid rows
pitch_physics$shape_cluster[valid_idx] <- predict(
  pitch_gmm,
  newdata = pitch_matrix_full[valid_idx, ]
)$classification

pitch_physics$shape_cluster <- factor(pitch_physics$shape_cluster)

# --------------------------------------------------
# CLUSTER CENTERS (CRITICAL FOR INTERPRETATION)
# --------------------------------------------------

pitch_cluster_centers <- pitch_physics %>%
  filter(!is.na(shape_cluster)) %>%
  group_by(shape_cluster) %>%
  summarise(
    velo = mean(velo_raw),
    spin = mean(spin_raw),
    ivb  = mean(ivb_raw),
    hb   = mean(hb_raw),
    n = n(),
    .groups = "drop"
  )

print(pitch_cluster_centers)



# --------------------------------------------------
# QUICK CLUSTER VISUALIZATION (HB vs IVB)
# --------------------------------------------------

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


# --------------------------------------------------
# BUILD PITCHER MIXTURE (PITCH SHAPE USAGE)
# --------------------------------------------------

pitch_mix <- pitch_physics %>%
  filter(!is.na(shape_cluster)) %>%
  group_by(pitcher_id, year, shape_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pitcher_id, year) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

# Convert to wide format
pitch_mix_wide <- pitch_mix %>%
  select(pitcher_id, year, shape_cluster, pct) %>%
  pivot_wider(
    names_from = shape_cluster,
    values_from = pct,
    values_fill = 0,
    names_prefix = "cluster_"
  )

# --------------------------------------------------
# ADD ENTROPY (ARSENAL DIVERSITY)
# --------------------------------------------------

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

# --------------------------------------------------
# PREP MATRIX FOR PITCHER CLUSTERING
# --------------------------------------------------

pitch_matrix <- pitch_features %>%
  select(starts_with("cluster_"), entropy) %>%
  mutate(across(where(is.numeric),
                ~ifelse(is.infinite(.), NA, .))) %>%
  drop_na()

# --------------------------------------------------
# CLUSTER PITCHERS (ARCHETYPES) - NO PCA
# --------------------------------------------------

set.seed(42)

gmm_pitchers <- Mclust(
  pitch_matrix,
  G = 3:8
)

pitcher_clusters <- pitch_features %>%
  select(pitcher_id, year) %>%
  distinct() %>%
  mutate(cluster = factor(gmm_pitchers$classification))

# --------------------------------------------------
# INSPECT PITCHER CLUSTERS (VERY IMPORTANT)
# --------------------------------------------------

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

# --------------------------------------------------
# BUILD PA-LEVEL DATA
# --------------------------------------------------

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

# --------------------------------------------------
# JOIN PITCHER CLUSTERS TO PA DATA
# --------------------------------------------------
pitcher_clusters <- pitcher_clusters %>%
  mutate(pitcher_id = as.character(pitcher_id))

pa_data_pitcher <- pa_data %>%
  mutate(pitcher_id = as.character(pitcher_id)) %>%
  inner_join(pitcher_clusters, by = c("pitcher_id","year"))

# --------------------------------------------------
# COORS INTERACTION MODEL (PITCHER ARCHETYPES)
# --------------------------------------------------

pitcher_cluster_model <- lmer(
  runs ~ coors * cluster + rockies + factor(year) +
    (1 | pitcher_id) +
    (1 | batter_id),
  data = pa_data_pitcher,
  REML = FALSE
)

summary(pitcher_cluster_model)

# --------------------------------------------------
# BUILD PITCH SHAPE USAGE FEATURES (FOR TRAIT MODEL)
# --------------------------------------------------

pitch_usage <- pitch_mix %>%
  select(pitcher_id, year, shape_cluster, pct) %>%
  pivot_wider(
    names_from = shape_cluster,
    values_from = pct,
    values_fill = 0,
    names_prefix = "pitch_"
  )

# Ensure type consistency
pitch_usage <- pitch_usage %>%
  mutate(pitcher_id = as.character(pitcher_id))

# --------------------------------------------------
# JOIN PITCH USAGE TO PA DATA
# --------------------------------------------------

pa_data_usage <- pa_data %>%
  mutate(pitcher_id = as.character(pitcher_id)) %>%
  inner_join(pitch_usage, by = c("pitcher_id","year"))

# --------------------------------------------------
# OPTIONAL: REMOVE ONE PITCH TO AVOID PERFECT COLLINEARITY
# (because all pitch_% sum to 1)
# --------------------------------------------------

# Drop pitch_1 (gyro slider) as baseline
pa_data_usage <- pa_data_usage %>%
  select(-pitch_1)

# --------------------------------------------------
# PITCH SHAPE TRAIT INTERACTION MODEL
# --------------------------------------------------

pitch_trait_model <- lmer(
  runs ~ coors * (
    pitch_2 + pitch_3 + pitch_4 +
      pitch_5 + pitch_6 + pitch_7 + pitch_8
  ) + 
    rockies +
    factor(year) +
    (1 | pitcher_id) +
    (1 | batter_id),
  data = pa_data_usage,
  REML = FALSE
)

summary(pitch_trait_model)

# --------------------------------------------------
# ADD ENTROPY TO PA DATA
# --------------------------------------------------

pitch_entropy <- pitch_entropy %>%
  mutate(pitcher_id = as.character(pitcher_id))

pa_data_entropy <- pa_data %>%
  mutate(pitcher_id = as.character(pitcher_id)) %>%
  inner_join(pitch_entropy, by = c("pitcher_id","year"))

# --------------------------------------------------
# ENTROPY INTERACTION MODEL
# --------------------------------------------------

entropy_model <- lmer(
  runs ~ coors * entropy +
    rockies +
    factor(year) +
    (1 | pitcher_id) +
    (1 | batter_id),
  data = pa_data_entropy,
  REML = FALSE
)

summary(entropy_model)

pitch_bucket_mix <- pitch_physics %>%
  group_by(pitcher_id, year, pitch_bucket) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pitcher_id, year) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

pitch_bucket_wide <- pitch_bucket_mix %>%
  select(pitcher_id, year, pitch_bucket, pct) %>%
  pivot_wider(
    names_from = pitch_bucket,
    values_from = pct,
    values_fill = 0,
    names_prefix = "pitch_"
  )

pitch_bucket_wide <- pitch_bucket_wide %>%
  mutate(pitcher_id = as.character(pitcher_id))

pa_data_bucket <- pa_data %>%
  mutate(pitcher_id = as.character(pitcher_id)) %>%
  inner_join(pitch_bucket_wide, by = c("pitcher_id","year"))

pa_data_bucket <- pa_data_bucket %>%
  select(-pitch_FF)  # or whatever baseline you want

bucket_model <- lmer(
  runs ~ coors * (
    pitch_SI + pitch_SL + pitch_CURVE +
      pitch_SWEEPER + pitch_CH + pitch_SPLIT + pitch_FC
  ) +
    rockies +
    factor(year) +
    (1 | pitcher_id) +
    (1 | batter_id),
  data = pa_data_bucket,
  REML = FALSE
)

summary(bucket_model)
