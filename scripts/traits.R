library(tidyverse)
library(data.table)
library(lubridate)
library(parallel)
library(lme4)
library(mclust)
library(sabRmetrics)
library(circular)
library(MuMIn)
library(gt)
library(baseballr)
library(stringr)
library(here)

data_path <- here("data", "savant_data")

years <- 2015:2025

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
    fread(f, select = cols_keep, showProgress = TRUE)
  }),
  use.names = TRUE, fill = TRUE
)

statcast <- statcast %>%
  mutate(
    game_date = as.Date(game_date),
    year = as.integer(year),
    coors = ifelse(home_team == "COL", 1, 0),
    rockies = ifelse(
      (home_team == "COL" & inning_topbot == "Bot") |
        (away_team == "COL" & inning_topbot == "Top"), 1, 0)
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
    rockies = first(rockies),
    runs = sum(delta_run_exp, na.rm = TRUE),
    .groups = "drop"
  )

pa_data <- pa_data %>%
  mutate(batter_id = as.character(batter_id))

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

spray_angle_savant <- function(hx, hy) {
  raw <- atan((hx - 125.42) / (198.27 - hy)) * 180/pi * 0.75
  (raw + 360) %% 360
}

ang_dist <- function(a, b) abs((a - b + 180) %% 360 - 180)

batted_balls <- statcast_noncoors %>%
  filter(!is.na(launch_angle),
         !is.na(launch_speed),
         !is.na(hit_coord_x),
         !is.na(hit_coord_y)) %>%
  baseballr::code_barrel() %>% 
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
    
    LA_consistency = sd_LA / mean_EV,
    
    spray_var = var.circular(circular(spray_deg, units="degrees")) / (2*pi),
    
    sweet_spot_rate = mean(sweet_spot, na.rm = TRUE),
    hard_hit_rate = mean(hard_hit, na.rm = TRUE),
    
    mean_EV_air = mean(launch_speed[air_ball], na.rm = TRUE),
    p90_EV_air = ifelse(sum(air_ball) > 0,
                        quantile(launch_speed[air_ball], 0.9, na.rm = TRUE),
                        NA_real_),
    
    gb_rate = mean(la_bucket == "ground", na.rm = TRUE),
    ld_rate = mean(la_bucket == "line", na.rm = TRUE),
    fb_rate = mean(la_bucket == "fly", na.rm = TRUE),
    
    pull_rate = mean(pull, na.rm = TRUE),
    oppo_rate = mean(oppo, na.rm = TRUE),
    pull_air_rate = ifelse(sum(air_ball) > 0, sum(pull & air_ball) / sum(air_ball), NA_real_),
    
    oppo_air_rate = ifelse(sum(air_ball) > 0, sum(oppo & air_ball) / sum(air_ball), NA_real_),
    
    pull_hard_hit_rate = ifelse(sum(pull) > 0, sum(pull & hard_hit) / sum(pull), NA_real_),
    
    oppo_hard_hit_rate = ifelse(sum(oppo) > 0, sum(oppo & hard_hit) / sum(oppo), NA_real_),
    
    barrel_rate = mean(barrel, na.rm = TRUE),
    
    .groups = "drop"
  )

hitter_season_traits <- batted_traits %>%
  inner_join(discipline_traits,
             by = c("batter_id","year")) %>%
  filter(n_bbe >= 150)

sprint_path <- here("data", "sprint_speeds")

sprint_files <- list.files(
  path = sprint_path,
  pattern = "sprint_speed_.*\\.csv",
  full.names = TRUE
)

if (length(sprint_files) == 0) {
  stop("No sprint speed files found in data/sprint_speeds")
}

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

trait_data_raw <- hitter_season_traits %>%
  select(
    batter_id,
    year,
    mean_LA,
    p90_EV,
    gb_rate,
    ld_rate,
    fb_rate,
    mean_EV_air,
    pull_air_rate,
    oppo_air_rate,
    barrel_rate,
    spray_var,
    hard_hit_rate,
    sweet_spot_rate
  ) %>%
  drop_na() %>%
  rename(
    mean_LA_raw = mean_LA,
    p90_EV_raw = p90_EV,
    gb_rate_raw = gb_rate,
    ld_rate_raw = ld_rate,
    fb_rate_raw = fb_rate,
    mean_EV_air_raw = mean_EV_air,
    pull_air_rate_raw = pull_air_rate,
    oppo_air_rate_raw = oppo_air_rate,
    barrel_rate_raw = barrel_rate,
    spray_var_raw = spray_var,
    hard_hit_rate_raw = hard_hit_rate,
    sweet_spot_rate_raw = sweet_spot_rate
  )

trait_data <- trait_data_raw %>%
  mutate(
    mean_LA = as.numeric(scale(mean_LA_raw)),
    p90_EV = as.numeric(scale(p90_EV_raw)),
    gb_rate = as.numeric(scale(gb_rate_raw)),
    ld_rate = as.numeric(scale(ld_rate_raw)),
    fb_rate = as.numeric(scale(fb_rate_raw)),
    mean_EV_air = as.numeric(scale(mean_EV_air_raw)),
    pull_air_rate = as.numeric(scale(pull_air_rate_raw)),
    oppo_air_rate = as.numeric(scale(oppo_air_rate_raw)),
    barrel_rate = as.numeric(scale(barrel_rate_raw)),
    spray_var = as.numeric(scale(spray_var_raw)),
    hard_hit_rate = as.numeric(scale(hard_hit_rate_raw)),
    sweet_spot_rate = as.numeric(scale(sweet_spot_rate_raw))
  ) %>%
  select(
    batter_id,
    year,
    mean_LA,
    p90_EV,
    gb_rate,
    ld_rate,
    fb_rate,
    mean_EV_air,
    pull_air_rate,
    oppo_air_rate,
    barrel_rate,
    spray_var,
    hard_hit_rate,
    sweet_spot_rate
  )

pa_data_traits <- pa_data %>%
  inner_join(trait_data, by = c("batter_id","year"))

continuous_model <- lmer(
  runs ~ coors *
    (p90_EV +
       gb_rate +
       fb_rate +
       mean_EV_air +
       barrel_rate +
       pull_air_rate) +
    rockies +
    factor(year) +
    (1 | batter_id) +
    (1 | pitcher_id),
  data = pa_data_traits,
  REML = FALSE
)

summary(continuous_model)

r.squaredGLMM(continuous_model)

baseline_model <- lmer(
  runs ~ coors +
    rockies +
    factor(year) +
    (1 | batter_id) +
    (1 | pitcher_id),
  data = pa_data_traits,
  REML = FALSE
)

anova_res <- anova(baseline_model, continuous_model)
anova_res

model_comparison <- tibble::tibble(
  Model = c("Baseline", "Trait Model"),
  AIC = anova_res$AIC,
  ChiSq = c(NA, anova_res$Chisq[2]),
  p_value = c(NA, anova_res$`Pr(>Chisq)`[2])
)

model_comparison <- model_comparison %>%
  mutate(
    p_value = ifelse(
      is.na(p_value),
      NA,
      format.pval(p_value, digits = 3, eps = 1e-16)
    )
  )

model_gt <- model_comparison %>%
  gt() %>%
  cols_label(
    Model = "Model",
    AIC = "AIC",
    ChiSq = "χ²",
    p_value = "p-value"
  ) %>%
  fmt_number(
    columns = c(AIC, ChiSq),
    decimals = 0
  ) %>%
  cols_align(
    align = "center",
    -Model
  ) %>%
  tab_header(
    title = "Model Comparison",
    subtitle = "Trait Model vs Baseline"
  )

model_gt

gtsave(model_gt, "model_comparison_table_clean.png")

fe <- fixef(continuous_model)

player_lookup <- statcast %>%
  select(batter_id, batter_name) %>%
  distinct() %>%
  mutate(batter_id = as.character(batter_id))

coef_df <- data.frame(
  feature = names(fe),
  value = fe
) %>%
  filter(str_detect(feature, "coors:")) %>%
  mutate(
    feature_clean = case_when(
      feature == "coors:p90_EV" ~ "90th Percentile EV",
      feature == "coors:fb_rate" ~ "Fly Ball Rate",
      feature == "coors:mean_EV_air" ~ "Avg EV on LDs and FBs",
      feature == "coors:barrel_rate" ~ "Barrel Rate",
      feature == "coors:pull_air_rate" ~ "Pull Air Rate",
      feature == "coors:gb_rate" ~ "Ground Ball Rate",
      TRUE ~ feature
    )
  )

feature_importance <- ggplot(
  coef_df,
  aes(x = reorder(feature_clean, value), y = value)
) +
  geom_col(width = 0.7, fill = "#4C78A8") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_text(
    aes(label = round(value, 3)),
    hjust = ifelse(value > 0, -0.1, 1.1)
    size = 4
  ) +
  coord_flip() +
  labs(
    title = "Trait Effects on Coors Advantage",
    subtitle = "Effect of 1 SD Increase in Trait on Run Value at Coors",
    x = "",
    y = "Change in Run Value per PA"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12)
  )

feature_importance
ggsave("feature_importance_clean.png", feature_importance, width = 14, height = 5)
