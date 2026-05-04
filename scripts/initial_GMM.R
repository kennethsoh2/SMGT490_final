library(tidyverse)
library(data.table)
library(lubridate)
library(parallel)
library(lme4)
library(mclust)
library(sabRmetrics)

setwd("~/R/SMGT 490")

years <- 2015:2025

# for (y in years) {
#   cl <- makeCluster(detectCores())
#   assign(
#     paste0("data_", y),
#     download_baseballsavant(
#       start_date = paste0(y, "-01-01"),
#       end_date   = paste0(y, "-12-31"),
#       cl = cl
#     )
#   )
#   stopCluster(cl)
# }

cols_keep <- c(
  "game_id","game_date","year","game_type",
  "home_team","away_team",
  "inning","inning_topbot","outs",
  "at_bat_number","event_index",
  "batter_id","batter_name","bat_side",
  "pitcher_id","pitch_hand",
  "launch_speed","launch_angle",
  "hit_coord_x","hit_coord_y",
  "plate_x","plate_z","zone",
  "description","balls","strikes",
  "events","woba_value","woba_denom",
  "delta_run_exp","delta_home_win_exp"
)

data_path <- "savant_data"

files <- list.files(
  path = data_path,
  pattern = "data_.*\\.csv",
  full.names = TRUE
)

statcast_list <- lapply(files, function(f) {
  fread(f, select = cols_keep, showProgress = TRUE)
})

statcast <- rbindlist(statcast_list, use.names = TRUE, fill = TRUE)
rm(statcast_list)

statcast <- statcast %>%
  mutate(
    game_date = as.Date(game_date),
    year = as.integer(year),
    coors = ifelse(home_team == "COL", 1, 0)
  ) %>%
  filter(game_type == "R")

pa_data <- statcast %>%
  filter(!is.na(woba_value), woba_denom > 0) %>%
  group_by(game_id, at_bat_number) %>%
  summarise(
    batter_id = first(batter_id),
    pitcher_id = first(pitcher_id),
    year = first(year),
    coors = first(coors),
    runs = sum(delta_run_exp, na.rm = TRUE),
    .groups = "drop"
  )

initial_model <- lmer(
  runs ~ coors + factor(year) + (1 | batter_id) + (1 | pitcher_id) + (1 + coors | batter_id:year),
  data = pa_data
)

ranef_slopes <- ranef(initial_model)$`batter_id:year`

player_effects <- ranef_slopes %>%
  rownames_to_column("batter_season") %>%
  separate(batter_season, into = c("batter_id","year"), sep=":") %>%
  rename(coors_re = coors) %>%
  mutate(
    year = as.integer(year),
    coors_effect = fixef(initial_model)["coors"] + coors_re
  ) %>%
  select(batter_id, year, coors_effect)

spray_angle_savant <- function(hx, hy) {
  raw <- atan((hx - 125.42) / (198.27 - hy)) * 180/pi * 0.75
  (raw + 360) %% 360
}

ang_dist <- function(a, b) abs((a - b + 180) %% 360 - 180)

statcast_noncoors <- statcast %>%
  filter(coors == 0)

pa_discipline <- statcast_noncoors %>%
  group_by(game_id, at_bat_number) %>%
  summarise(
    batter_id = first(batter_id),
    year = first(year),
    is_strikeout = any(events == "strikeout", na.rm = TRUE),
    is_walk = any(events %in% c("walk","intent_walk"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(batter_id, year) %>%
  summarise(
    n_PA = n(),
    K_rate = sum(is_strikeout) / n_PA,
    BB_rate = sum(is_walk) / n_PA,
    .groups = "drop"
  )

pitch_discipline <- statcast_noncoors %>%
  mutate(
    in_zone = zone %in% 1:9,
    
    swing = description %in% c(
      "swinging_strike",
      "swinging_strike_blocked",
      "foul",
      "hit_into_play"
    ),
    
    contact = description %in% c(
      "foul",
      "hit_into_play"
    ),
    
    whiff = description %in% c(
      "swinging_strike",
      "swinging_strike_blocked"
    )
  ) %>%
  group_by(batter_id, year) %>%
  summarise(
    total_pitches = n(),
    
    swinging_strike_rate = sum(whiff) / total_pitches,
    
    zone_swings = sum(in_zone & swing),
    zone_contacts = sum(in_zone & contact),
    zone_whiffs = sum(in_zone & whiff),
    
    zone_contact_rate = ifelse(zone_swings > 0,
                               zone_contacts / zone_swings,
                               NA_real_),
    
    zone_whiff_rate = ifelse(zone_swings > 0,
                             zone_whiffs / zone_swings,
                             NA_real_),
    
    .groups = "drop"
  )

discipline_traits <- pitch_discipline %>%
  inner_join(pa_discipline,
             by = c("batter_id","year"))

batted_balls <- statcast_noncoors %>%
  filter(!is.na(launch_angle),
         !is.na(launch_speed),
         !is.na(hit_coord_x),
         !is.na(hit_coord_y)) %>%
  mutate(
    spray_deg = spray_angle_savant(hit_coord_x, hit_coord_y),
    dist_to_lf = ang_dist(spray_deg, 315),
    dist_to_rf = ang_dist(spray_deg, 45),
    pull = case_when(
      bat_side == "R" ~ dist_to_lf <= 20,
      bat_side == "L" ~ dist_to_rf <= 20,
      TRUE ~ FALSE
    ),
    oppo = case_when(
      bat_side == "R" ~ dist_to_rf <= 20,
      bat_side == "L" ~ dist_to_lf <= 20,
      TRUE ~ FALSE
    ),
    la_bucket = case_when(
      launch_angle < 10 ~ "ground",
      launch_angle >= 10 & launch_angle < 25 ~ "line",
      launch_angle >= 25 & launch_angle <= 50 ~ "fly",
      launch_angle > 50 ~ "popup",
      TRUE ~ NA_character_
    ),
    sweet_spot = launch_angle >= 8 & launch_angle <= 32,
    air_ball = launch_angle >= 10 & launch_angle <= 50,
    hard_hit = launch_speed >= 95
  )

batted_traits <- batted_balls %>%
  group_by(batter_id, year) %>%
  summarise(
    n_bbe = n(),
    
    mean_LA = mean(launch_angle, na.rm = TRUE),
    sd_LA = sd(launch_angle, na.rm = TRUE),
    IQR_LA = IQR(launch_angle, na.rm = TRUE),
    
    mean_EV = mean(launch_speed, na.rm = TRUE),
    p90_EV = quantile(launch_speed, 0.9, na.rm = TRUE),
    sd_EV = sd(launch_speed, na.rm = TRUE),
    IQR_EV = IQR(launch_speed, na.rm = TRUE),
    
    sweet_spot_rate = mean(sweet_spot, na.rm = TRUE),
    
    mean_EV_air = mean(launch_speed[air_ball], na.rm = TRUE),
    p90_EV_air = quantile(launch_speed[air_ball], 0.9, na.rm = TRUE),
    
    gb_rate = mean(la_bucket == "ground", na.rm = TRUE),
    ld_rate = mean(la_bucket == "line", na.rm = TRUE),
    fb_rate = mean(la_bucket == "fly", na.rm = TRUE),
    
    pull_rate = mean(pull, na.rm = TRUE),
    oppo_rate = mean(oppo, na.rm = TRUE),
    pull_air_rate = mean(pull & air_ball, na.rm = TRUE),
    oppo_air_rate = mean(oppo & air_ball, na.rm = TRUE),
    pull_hard_hit_rate = mean(pull & hard_hit, na.rm = TRUE),
    oppo_hard_hit_rate = mean(oppo & hard_hit, na.rm = TRUE),
    
    .groups = "drop"
  )

hitter_season_traits <- batted_traits %>%
  inner_join(discipline_traits,
             by = c("batter_id","year")) %>%
  filter(n_bbe >= 150)

sprint_path <- "sprint_speeds"

sprint_files <- list.files(
  path = sprint_path,
  pattern = "sprint_speed_.*\\.csv",
  full.names = TRUE
)

years_extracted <- as.integer(
  str_extract(basename(sprint_files), "\\d{4}")
)

sprint_list <- Map(function(f, y) {
  
  dt <- fread(
    f,
    select = c(
      "last_name, first_name",
      "player_id",
      "team",
      "sprint_speed"
    ),
    showProgress = TRUE
  )
  
  dt[, year := y]
  dt
  
}, sprint_files, years_extracted)

sprint_df <- rbindlist(sprint_list, use.names = TRUE, fill = TRUE)

rm(sprint_list)

sprint_df <- sprint_df %>%
  rename(
    batter_id = player_id,
    player_name = `last_name, first_name`
  ) %>%
  mutate(
    batter_id = as.character(batter_id),
    year = as.integer(year)
  )

hitter_season_traits <- hitter_season_traits %>%
  mutate(
    batter_id = as.character(batter_id)
  )

hitter_season_traits <- hitter_season_traits %>%
  left_join(
    sprint_df,
    by = c("batter_id", "year")
  )

#PCA
trait_data <- hitter_season_traits %>%
  select(batter_id, year, where(is.numeric)) %>%
  select(-n_bbe) %>%
  drop_na()

trait_matrix <- trait_data %>%
  select(-batter_id, -year)

trait_scaled <- scale(trait_matrix)
pca_model <- prcomp(trait_scaled)

var_exp <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)
num_pc <- which(var_exp >= 0.80)[1]

pca_scores <- as.data.frame(pca_model$x[,1:num_pc])

#GMM
gmm_model <- Mclust(pca_scores)

posterior_df <- as.data.frame(gmm_model$z)
colnames(posterior_df) <- paste0("Cluster",1:ncol(posterior_df))

posterior_df$batter_id <- trait_data$batter_id
posterior_df$year <- trait_data$year

#Regress Coors effect on cluster weights
cluster_df <- posterior_df %>%
  inner_join(player_effects,
             by = c("batter_id","year"))

cluster_vars <- colnames(cluster_df)[
  grepl("Cluster", colnames(cluster_df))
]

cluster_vars <- cluster_vars[-length(cluster_vars)]

formula_str <- paste(
  "coors_effect ~",
  paste(cluster_vars, collapse = " + ")
)

cluster_model <- lm(as.formula(formula_str),
                    data = cluster_df)

summary(cluster_model)

cluster_summary <- posterior_df %>%
  mutate(cluster = gmm_model$classification) %>%
  left_join(trait_data, by = c("batter_id","year")) %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean),
            n = n(),
            .groups = "drop")

write.csv(cluster_summary, "cluster_summary_2.csv")

