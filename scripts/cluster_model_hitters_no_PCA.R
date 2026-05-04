library(tidyverse)
library(gt)
library(data.table)
library(lubridate)
library(parallel)
library(lme4)
library(mclust)
library(sabRmetrics)
library(baseballr)
library(circular)
library(broom.mixed)

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
  "delta_run_exp","delta_home_win_exp",
  "bat_speed", "swing_length", "swing_path_tilt", "intercept_ball_minus_batter_pos_y_inches"
)

data_path <- "savant_data"

files <- list.files(
  path = data_path,
  pattern = "data_.*\\.csv",
  full.names = TRUE
)

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
    
    swinging_strike_rate = ifelse(total_pitches > 0,
                                  sum(whiff) / total_pitches,
                                  NA_real_),
    
    zone_swings = sum(in_zone & swing),
    zone_contacts = sum(in_zone & contact),
    zone_whiffs = sum(in_zone & whiff),
    
    zone_swing_rate = ifelse(sum(in_zone) > 0,
                             sum(in_zone & swing) / sum(in_zone),
                             NA_real_),
    
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
    p90_EV_air = quantile(launch_speed[air_ball], 0.9, na.rm = TRUE),
    
    gb_rate = mean(la_bucket == "ground", na.rm = TRUE),
    ld_rate = mean(la_bucket == "line", na.rm = TRUE),
    fb_rate = mean(la_bucket == "fly", na.rm = TRUE),
    
    pull_rate = mean(pull, na.rm = TRUE),
    oppo_rate = mean(oppo, na.rm = TRUE),
    pull_air_rate = ifelse(sum(air_ball) > 0,
                           sum(pull & air_ball) / sum(air_ball),
                           NA_real_),
    
    oppo_air_rate = ifelse(sum(air_ball) > 0,
                           sum(oppo & air_ball) / sum(air_ball),
                           NA_real_),
    
    pull_hard_hit_rate = ifelse(sum(pull) > 0,
                                sum(pull & hard_hit) / sum(pull),
                                NA_real_),
    
    oppo_hard_hit_rate = ifelse(sum(oppo) > 0,
                                sum(oppo & hard_hit) / sum(oppo),
                                NA_real_),
    
    barrel_rate = mean(barrel, na.rm = TRUE),
    
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

# --------------------------------------------------
# PREP TRAIT DATA
# --------------------------------------------------

trait_data <- hitter_season_traits %>%
  select(
    -n_bbe,
    -n_PA,
    -total_pitches,
    -zone_swings,
    -zone_contacts,
    -zone_whiffs
  )

trait_data_clean <- trait_data %>%
  mutate(across(where(is.numeric),
                ~ifelse(is.infinite(.), NA, .))) %>%
  drop_na()

trait_matrix <- trait_data_clean %>% 
  select(
    mean_LA, sd_LA, LA_consistency,
    p90_EV, barrel_rate,
    gb_rate, fb_rate,
    pull_rate, pull_air_rate, spray_var,
    zone_swing_rate, zone_contact_rate, swinging_strike_rate,
    K_rate, BB_rate,
    sprint_speed
  )

# --------------------------------------------------
# SCALE + GMM (NO PCA)
# --------------------------------------------------

set.seed(42)

trait_matrix_scaled <- scale(trait_matrix)

gmm_model_2 <- Mclust(trait_matrix_scaled, G = 3:6)

# Check cluster sizes
table(gmm_model_2$classification)

# BIC plot (correct model)
plot(gmm_model_2, what = "BIC")

# --------------------------------------------------
# ASSIGN CLUSTERS
# --------------------------------------------------

cluster_assignments <- data.frame(
  batter_id = trait_data_clean$batter_id,
  year = trait_data_clean$year,
  cluster = factor(gmm_model_2$classification)
)

# --------------------------------------------------
# JOIN TO PA DATA
# --------------------------------------------------

pa_data_clustered <- pa_data %>%
  inner_join(cluster_assignments, by = c("batter_id","year"))

# --------------------------------------------------
# MIXED EFFECTS MODEL
# --------------------------------------------------

cluster_model <- lmer(
  runs ~ coors * cluster + rockies + factor(year) +
    (1 | batter_id) +
    (1 | pitcher_id),
  data = pa_data_clustered,
  REML = FALSE
)

summary(cluster_model)

# --------------------------------------------------
# CLEAN MODEL OUTPUT TABLE
# --------------------------------------------------

model_results <- broom.mixed::tidy(cluster_model, effects = "fixed")

model_results_clean <- model_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_when(
      term == "coors" ~ "Coors Field",
      grepl("^cluster", term) ~ gsub("cluster", "Cluster ", term),
      grepl("^coors:cluster", term) ~ gsub("coors:cluster", "Coors × Cluster ", term),
      TRUE ~ term
    )
  ) %>%
  filter(grepl("Coors|Cluster", term))

model_results_clean <- model_results_clean %>%
  mutate(
    estimate = round(estimate, 4),
    std.error = round(std.error, 4),
    t.value = round(statistic, 2)
  ) %>% 
  select(term, estimate, std.error, t.value)

model_results_table <- model_results_clean %>%
  gt() %>%
  cols_label(
    term = "Variable",
    estimate = "Coefficient",
    std.error = "Std. Error",
    t.value = "t-value"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_header(
    title = "Mixed Effects Model: Coors Field × Hitter Archetype"
  )

model_results_table
gtsave(model_results_table, "cluster_model_results_v2.png")

# --------------------------------------------------
# CLUSTER SUMMARY TABLE
# --------------------------------------------------

cluster_summary <- trait_data_clean %>%
  left_join(cluster_assignments, by = c("batter_id","year")) %>%
  left_join(
    statcast %>%
      distinct(batter_id = as.character(batter_id), batter_name),
    by = "batter_id"
  ) %>%
  group_by(cluster) %>%
  summarise(
    across(where(is.numeric) & !year, mean),
    example_players = paste(head(unique(batter_name), 5), collapse = ", "),
    n = n(),
    .groups = "drop"
  )

cluster_summary_final <- cluster_summary %>% 
  select(cluster, mean_EV, p90_EV, pull_rate, pull_air_rate, barrel_rate, 
         zone_swing_rate, zone_contact_rate, K_rate, BB_rate, sprint_speed) %>%
  gt() %>%
  cols_label(
    cluster = "Cluster",
    mean_EV = "Avg EV",
    p90_EV = "90th Percentile EV",
    pull_rate = "Pull%",
    pull_air_rate = "Air Pull%",
    barrel_rate = "Barrel%",
    zone_swing_rate = "Zone Swing%",
    zone_contact_rate = "Zone Contact%",
    K_rate = "K%",
    BB_rate = "BB%",
    sprint_speed = "Sprint Speed"
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  fmt_number(
    columns = c(mean_EV, p90_EV, sprint_speed),
    decimals = 2
  ) %>% 
  fmt_percent(
    columns = c(pull_rate, pull_air_rate, barrel_rate,
                zone_swing_rate, zone_contact_rate,
                K_rate, BB_rate),
    decimals = 2
  ) %>%
  data_color(
    columns = c(pull_rate, pull_air_rate, barrel_rate,
                zone_swing_rate, zone_contact_rate,
                BB_rate),
    colors = scales::col_numeric(
      palette = c("blue", "white", "red"),
      domain = NULL
    )
  ) %>%
  data_color(
    columns = c(mean_EV, p90_EV, sprint_speed),
    colors = scales::col_numeric(
      palette = c("blue", "white", "red"),
      domain = NULL
    )
  ) %>%
  data_color(
    columns = K_rate,
    colors = scales::col_numeric(
      palette = c("red", "white", "blue"),
      domain = NULL
    )
  ) %>%
  tab_header(
    title = "Hitter Archetype Clusters"
  )

cluster_summary_final
gtsave(cluster_summary_final, "cluster_summary_final_no_pca.png")

# --------------------------------------------------
# MODEL COMPARISON
# --------------------------------------------------

model_no_interaction <- lmer(
  runs ~ coors + cluster + rockies + factor(year) +
    (1 | batter_id) +
    (1 | pitcher_id),
  data = pa_data_clustered,
  REML = FALSE
)

anova(model_no_interaction, cluster_model)

anova_results <- anova(model_no_interaction, cluster_model)

anova_df <- as.data.frame(anova_results) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(
    Model = c("Baseline Model", "Trait Interaction Model"),
    Parameters = Df,
    AIC = round(AIC, 0),
    BIC = round(BIC, 0),
    `Log Likelihood` = round(logLik, 0),
    `Chi-Square` = round(Chisq, 0),
    df = Df,
    `p-value` = case_when(
      is.na(`Pr(>Chisq)`) ~ NA_character_,
      `Pr(>Chisq)` < 2.2e-16 ~ "< 2.2e-16",
      TRUE ~ formatC(`Pr(>Chisq)`, format = "e", digits = 2)
    )
  ) %>%
  select(
    Model,
    Parameters,
    AIC,
    BIC,
    `Log Likelihood`,
    `Chi-Square`,
    df,
    `p-value`
  )

anova_table <- anova_df %>%
  gt() %>%
  
  tab_header(
    title = "Model Comparison: Baseline vs Trait Interaction Model",
    subtitle = "Likelihood Ratio Test"
  ) %>%
  
  cols_label(
    Model = "Model",
    Parameters = "Parameters",
    AIC = "AIC",
    BIC = "BIC",
    `Log Likelihood` = "Log Likelihood",
    `Chi-Square` = "χ²",
    df = "df",
    `p-value` = "p-value"
  ) %>%
  
  cols_align(
    align = "center",
    columns = -Model
  ) %>%
  
  fmt_number(
    columns = c(AIC, BIC, `Log Likelihood`, `Chi-Square`),
    decimals = 0,
    use_seps = TRUE
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  
  tab_options(
    table.font.size = 14,
    data_row.padding = px(6)
  )

anova_table
gtsave(anova_table, "anova_model_comparison.png")

# --------------------------------------------------
# CLUSTER SIZE PLOT
# --------------------------------------------------

ggplot(
  as.data.frame(table(cluster_assignments$cluster)),
  aes(Var1, Freq, fill = Var1)
) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Cluster Sizes (Statcast Traits)",
    x = "Cluster",
    y = "Number of Batter-Seasons"
  ) +
  theme(legend.position = "none")

# --------------------------------------------------
# CLUSTER PROFILES
# --------------------------------------------------

cluster_profiles <- trait_data %>%
  left_join(cluster_assignments, by = c("batter_id","year")) %>%
  group_by(cluster) %>%
  summarise(
    mean_EV = mean(mean_EV),
    barrel_rate = mean(barrel_rate),
    hard_hit_rate = mean(hard_hit_rate),
    gb_rate = mean(gb_rate),
    pull_rate = mean(pull_rate),
    zone_contact_rate = mean(zone_contact_rate),
    swinging_strike_rate = mean(swinging_strike_rate),
    BB_rate = mean(BB_rate),
    K_rate = mean(K_rate),
    sprint_speed = mean(sprint_speed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-cluster)

ggplot(cluster_profiles,
       aes(name, value, fill = cluster)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(
    title = "Cluster Skill Profiles (Statcast Traits)",
    x = "Trait",
    y = "Average Value"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --------------------------------------------------
# COORS EFFECT BY CLUSTER
# --------------------------------------------------

coors_effects <- pa_data_clustered %>%
  group_by(cluster, coors) %>%
  summarise(mean_runs = mean(runs), .groups = "drop")

coors_effects_plot <- ggplot(
  coors_effects,
  aes(cluster, mean_runs, fill = factor(coors))
) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(
    title = "Run Value by Cluster (Statcast Traits)",
    x = "Cluster",
    y = "Average Run Value per PA",
    fill = "Coors"
  )

coors_effects_plot
ggsave("coors_effects_plot.png", coors_effects_plot)
