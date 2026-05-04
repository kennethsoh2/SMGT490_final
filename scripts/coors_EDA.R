library(tidyverse)
library(data.table)
library(lubridate)

#Getting 2015-2025 Statcast Data
if (!require("sabRmetrics", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("saberpowers/sabRmetrics")
}
library(sabRmetrics)

cluster <- parallel::makeCluster(parallel::detectCores())
data_2015 <- sabRmetrics::download_baseballsavant(
 start_date = "2015-01-01",
 end_date = "2015-12-31",
 cl = cluster
)
parallel::stopCluster(cluster)

write.csv(data_2015, "data_2015.csv")
#data_2015 <- fread("data_2015.csv")

#Variables of Interest
cols_keep <- c(
  "game_id", "game_date", "year", "game_type",
  "home_team", "away_team",
  "inning", "inning_topbot", "outs",
  "at_bat_number", "event_index",
  
  "batter_id", "batter_name", "bat_side",
  "pitcher_id", "pitch_hand",
  
  "pitch_type", "release_speed", "release_spin_rate",
  "pfx_x", "pfx_z", "extension", "spin_axis",
  
  "launch_speed", "launch_angle",
  "hit_coord_x", "hit_coord_y",
  "bb_type", "hit_distance_sc",
  
  "plate_x", "plate_z", "zone",
  "description", "balls", "strikes",
  
  "events", "woba_value", "woba_denom",
  "delta_run_exp", "delta_home_win_exp",
  "expected_woba",

  "bat_speed", "swing_length",
  "attack_angle", "attack_direction",
  "swing_path_tilt"
)

years <- 2015:2025

for (y in years) {
  obj_name <- paste0("data_", y)
  df <- get(obj_name)
  keep_cols_year <- intersect(cols_keep, names(df))
  df_reduced <- df[, keep_cols_year, with = FALSE]
  assign(obj_name, df_reduced)
  rm(df, df_reduced)
  gc()
}

data_list <- mget(paste0("data_", 2015:2025))
statcast <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

#Coors Indicator
statcast <- statcast %>%
  mutate(
    game_date = as.Date(game_date),
    year = as.integer(year),
    coors = ifelse(home_team == "COL", 1, 0),
    covid_season = ifelse(year == 2020, 1, 0),
    is_home_batter = ifelse(inning_topbot == "Bot", 1, 0)
  )
statcast <- statcast %>%
  filter(game_type == "R")

statcast <- statcast %>%
  mutate(
    batter_id = as.factor(batter_id),
    pitcher_id = as.factor(pitcher_id),
    pitch_type = as.factor(pitch_type),
    pitch_hand = as.factor(pitch_hand),
    bat_side = as.factor(bat_side),
    home_team = as.factor(home_team),
    away_team = as.factor(away_team)
  )

#Aggregate to PA-level
pa_data <- statcast %>%
  filter(!is.na(woba_value), woba_denom > 0) %>%
  group_by(game_id, at_bat_number) %>%
  summarise(
    year = first(year),
    coors = first(coors),
    runs = sum(delta_run_exp, na.rm = TRUE),
    woba_value = sum(woba_value, na.rm = TRUE),
    woba_denom = sum(woba_denom, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(woba = woba_value / woba_denom)

#Coors Effect: Run Environment
pa_data %>%
  group_by(coors) %>%
  summarise(
    avg_run_value = mean(runs, na.rm = TRUE),
    avg_woba = mean(woba, na.rm = TRUE),
    n = n()
  )

library(ggplot2)

pa_data %>%
  group_by(year, coors) %>%
  summarise(avg_run_value = mean(runs, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = year, y = avg_run_value, color = factor(coors))) +
  geom_line(size = 1) +
  labs(color = "Coors",
       title = "Average Run Value per PA by Year",
       y = "Avg Run Value")

#Within-Hitter Coors Effect Differences
hitter_pa <- pa_data %>%
  left_join(
    statcast %>%
      select(game_id, at_bat_number, batter_id) %>%
      distinct(),
    by = c("game_id", "at_bat_number")
  )

hitter_splits <- hitter_pa %>%
  group_by(batter_id, coors) %>%
  summarise(
    avg_run_value = mean(runs, na.rm = TRUE),
    avg_woba = mean(woba, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

hitter_wide <- hitter_splits %>%
  pivot_wider(
    names_from = coors,
    values_from = c(avg_run_value, avg_woba, n),
    names_prefix = "coors_"
  )

hitter_wide <- hitter_wide %>%
  mutate(
    run_diff = avg_run_value_coors_1 - avg_run_value_coors_0,
    woba_diff = avg_woba_coors_1 - avg_woba_coors_0
  )

hitter_wide_filtered <- hitter_wide %>%
  filter(
    n_coors_1 >= 50,
    n_coors_0 >= 200
  )

ggplot(hitter_wide_filtered, aes(x = run_diff)) +
  geom_histogram(bins = 40) +
  labs(title = "Distribution of Within-Hitter Coors Run Value Differences")

#Initial (simple) Mixed-Effects Model
library(lme4)

initial_model <- lmer(
  runs ~ coors + (1 + coors | batter_id),
  data = hitter_pa
)

summary(initial_model)
VarCorr(initial_model)

#Player-level Coors Effects
ranef_slopes <- ranef(initial_model)$batter_id
head(ranef_slopes)

player_effects <- ranef_slopes %>%
  tibble::rownames_to_column("batter_id") %>%
  rename(
    intercept_re = `(Intercept)`,
    coors_re = coors
  ) %>%
  mutate(coors_effect = fixef(initial_model)["coors"] + coors_re)

summary(player_effects$coors_effect)
sd(player_effects$coors_effect)

#Trait Construction
spray_angle_savant <- function(hx, hy) {
  raw <- atan((hx - 125.42) / (198.27 - hy)) * 180/pi * 0.75
  (raw + 360) %% 360  # wrap into 0–360
}

ang_dist <- function(a, b) abs((a - b + 180) %% 360 - 180)

statcast_noncoors <- statcast %>%
  filter(coors == 0)

batted_balls <- statcast_noncoors %>%
  filter(!is.na(launch_angle),
         !is.na(launch_speed),
         !is.na(hit_coord_x),
         !is.na(hit_coord_y)) %>%
  mutate(
    spray_deg_xy = spray_angle_savant(hit_coord_x, hit_coord_y)
  )

batted_balls <- batted_balls %>%
  mutate(
    dist_to_lf = ang_dist(spray_deg_xy, 315),
    dist_to_rf = ang_dist(spray_deg_xy, 45),
    pull = case_when(
      bat_side == "R" ~ dist_to_lf <= 20,
      bat_side == "L" ~ dist_to_rf <= 20,
      TRUE ~ FALSE
    ),
    oppo = case_when(
      bat_side == "R" ~ dist_to_rf <= 20,
      bat_side == "L" ~ dist_to_lf <= 20,
      TRUE ~ FALSE
    )
  )

hitter_traits <- batted_balls %>%
  group_by(batter_id) %>%
  summarise(
    n_bbe = n(),
    #Launch angle profile
    mean_LA = mean(launch_angle, na.rm = TRUE),
    sd_LA = sd(launch_angle, na.rm = TRUE),
    
    #Exit velocity profile
    mean_EV = mean(launch_speed, na.rm = TRUE),
    p90_EV = quantile(launch_speed, 0.9, na.rm = TRUE),
    
    #Batted ball distribution
    gb_rate = mean(bb_type == "ground_ball", na.rm = TRUE),
    fb_rate = mean(bb_type == "fly_ball", na.rm = TRUE),
    ld_rate = mean(bb_type == "line_drive", na.rm = TRUE),
    
    #Spray angle tendencies
    pull_rate = mean(pull, na.rm = TRUE),
    oppo_rate = mean(oppo, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  filter(n_bbe >= 200)

summary(hitter_traits)

#Cluster hitters by trait profiles
trait_slope_data <- hitter_traits %>%
  inner_join(player_effects, by = "batter_id")

trait_matrix <- hitter_traits %>%
  select(-batter_id, -n_bbe)

trait_scaled <- scale(trait_matrix)

set.seed(42)

wss <- sapply(2:8, function(k) {
  kmeans(trait_scaled, centers = k, nstart = 20)$tot.withinss
})

plot(2:8, wss, type = "b",
     xlab = "Number of clusters",
     ylab = "Within-cluster sum of squares")

set.seed(42)

km4 <- kmeans(trait_scaled, centers = 4, nstart = 50)

hitter_traits$cluster <- factor(km4$cluster)

cluster_summary <- hitter_traits %>%
  group_by(cluster) %>%
  summarise(
    across(mean_LA:oppo_rate, mean),
    n = n(),
    .groups = "drop"
  )

cluster_summary

#Clusters in Coors
cluster_coors <- hitter_traits %>%
  inner_join(player_effects, by = "batter_id") %>%
  group_by(cluster) %>%
  summarise(
    mean_coors_effect = mean(coors_effect),
    sd_coors_effect = sd(coors_effect),
    n = n(),
    .groups = "drop"
  )

cluster_coors
