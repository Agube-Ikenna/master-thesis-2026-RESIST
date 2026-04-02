# ==================================================================
# Leaf Decomposition in the Kinzig Catchment – Full R Workflow
# Author: Agube Ikenna Stephen
# Last updated: 04.07.2025
# Description: Mixed model analysis
# ==================================================================

# ------------------------------------------------------------
# 1. Setup: Load Libraries and Set Working Directory
# ------------------------------------------------------------

# ── Working directory (portable – no manual editing needed) ───────────────────
# Detects where this script is saved and points to the data/processed folder.
# Works on any machine as long as the folder structure is kept intact.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  script_dir <- getSrcDirectory(function() {})  # fallback for source()
}
data_dir <- normalizePath(file.path(script_dir, "..", "data", "processed"),
                          mustWork = FALSE)
setwd(data_dir)
message("Working directory: ", getwd())

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(janitor)
library(ggeffects)
library(patchwork)
library(broom.mixed)

# ------------------------------------------------------------
# 2. Load and Prepare All Field Data
# ------------------------------------------------------------

processed     <- read_excel("2025_field_processed_samples.xlsx")
initial_mass  <- read_excel("2025_field_initial_mass.xlsx")
transport     <- read_excel("2025_field_transport_controls.xlsx")
site_info     <- read_excel("2025_field_site_information_Kinzig.xlsx")
inverts       <- read_excel("2025_field_invertebrates_coarse_leaf_bags.xlsx")
phch_raw      <- read_excel("202505_phch_data_Kinzig.xlsx", sheet = 1) %>% clean_names()

# Process land use info
site_info <- site_info %>%
  dplyr::mutate(LandUse = case_when(
    tolower(trimws(rural)) == "x" ~ "rural",
    tolower(trimws(urban)) == "x" ~ "urban",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::select(stream, LandUse)

# Handling loss correction: Calculate mean initial AFDM from transport bags
transport <- transport %>% mutate(AFDM_0 = mass_60_g - mass_500_g)
mean_AFDM_0 <- mean(transport$AFDM_0, na.rm = TRUE)

# Recalculate AFDM and daily decomposition rate (k) normalized by degree days
kinzig_processed <- processed %>%
  dplyr::filter(catchment == "Kinzig") %>%
  dplyr::mutate(
    OMD_60_mg = as.numeric(OMD_60_mg),
    OMD_500_mg = as.numeric(OMD_500_mg),
    AFDM_final = OMD_60_mg - OMD_500_mg,
    AFDM_0 = mean_AFDM_0,             # Global mean correction
    deg_days = 14 * 12.5,             # 14 days × 12.5 °C (degree days)
    OMD_dday = ((AFDM_0 - AFDM_final) / AFDM_0) * 100 / deg_days
  )

# Merge land use info and filter complete cases
kinzig_labeled <- kinzig_processed %>%
  left_join(site_info, by = "stream") %>%
  filter(!is.na(LandUse), !is.na(OMD_dday))

# Reshape data: get fine (microbial), coarse (total), and macroinvertebrate (difference) decomposition rates
model_data <- kinzig_labeled %>%
  dplyr::select(stream, replicate, type, LandUse, OMD_dday) %>%
  pivot_wider(names_from = type, values_from = OMD_dday) %>%
  rename(k_micro = fine, k_total = coarse) %>%
  dplyr::mutate(k_macro = k_total - k_micro)

# ------------------------------------------------------------
# 3. Join Invertebrate Abundance Data
# ------------------------------------------------------------

inverts_clean <- inverts %>%
  dplyr::filter(type == "coarse") %>%
  dplyr::mutate(across(c(Gammarids, Trichoptera, Chironomids, Others), as.numeric)) %>%
  dplyr::mutate(total_abundance = rowSums(across(c(Gammarids, Trichoptera, Chironomids, Others)), na.rm = TRUE)) %>%
  dplyr::select(stream, replicate, total_abundance, Gammarids, Trichoptera, Chironomids, Others)

model_data <- model_data %>%
  dplyr::mutate(across(c(stream, replicate), str_trim)) %>%
  left_join(inverts_clean, by = c("stream", "replicate"))

# ------------------------------------------------------------
# 4. Merge Physico-Chemical Stream Parameters
# ------------------------------------------------------------

phch_clean <- phch_raw %>%
  dplyr::select(stream,
         ph = p_h,
         water_temp = water_temp,
         oxygen_dissolved = oxygen_dissolved,
         oxygen_saturation = oxygen_saturation,
         redox = redox_potential,
         conductivity = conduc_tivity,
         water_level = water_level_at_pole) %>%
  dplyr::filter(!is.na(stream)) %>%
  dplyr::mutate(stream = str_trim(stream))

model_data <- model_data %>%
  left_join(phch_clean, by = "stream")

# Check merged data structure
glimpse(model_data)

# ------------------------------------------------------------
# 5. Fit Linear Mixed-Effects Models
# ------------------------------------------------------------

# Hypothesis 1: Microbial decomposition ~ LandUse
model_micro <- lmer(k_micro ~ LandUse + (1 | stream), data = model_data)
summary(model_micro)
# INTERPRETATION: 
# The microbial decomposition rate (k_fine) showed no significant difference between 
# rural (M = -2.47 ± 0.09 %/degree-day) and urban (M = -2.43 ± 0.09 %/degree-day) 
# streams (t(8) = -0.32, p = 0.756). This suggests microbial communities were equally 
# efficient across land use types, possibly due to similar water chemistry conditions 
# supporting microbial activity.

# Hypothesis 2: Total decomposition ~ LandUse
model_total <- lmer(k_total ~ LandUse + (1 | stream), data = model_data)
summary(model_total)
# INTERPRETATION:
# Total decomposition rates (k_coarse) were 48.1% higher in rural streams (M = -0.81 ± 0.34) 
# compared to urban (M = -1.56 ± 0.35), though this difference was not statistically 
# significant (t(8) = 1.63, p = 0.142). The large effect size suggests potential biological 
# relevance despite limited statistical power from small sample sizes.

# Hypothesis 3: Macroinvertebrate-mediated decomposition ~ LandUse
model_macro <- lmer(k_macro ~ LandUse + (1 | stream), data = model_data)
summary(model_macro)
 #INTERPRETATION:
  # Macroinvertebrate-mediated decomposition (k_macro) showed a 91.4% increase in rural 
  # streams (1.66 ± 0.34) versus urban (0.87 ± 0.36), approaching significance (t(8) = 1.72, 
  # p = 0.128). This aligns with expectations of reduced shredder populations in urbanized 
  # streams due to habitat degradation and pollution sensitivity.
  
# Extended model including total macroinvertebrate abundance
model_macro_abund <- lmer(k_macro ~ LandUse + total_abundance + (1 | stream), data = model_data)
summary(model_macro_abund)
# Interpretation:
# Including abundance increased rural-urban difference (+0.90) and reduced p to 0.093, indicating biological relevance.

# Combined environmental covariates model
model_env <- lmer(k_macro ~ LandUse + water_temp + conductivity + oxygen_dissolved + (1 | stream), data = model_data)
summary(model_env)
# Interpretation:
# Oxygen significantly negatively affected k_macro (p=0.039), temperature marginal (p=0.066), conductivity no effect.


# ------------------------------------------------------------
# 6. Figures
# ------------------------------------------------------------

library(ggplot2)
library(ggeffects)
library(patchwork)

# Figure 1: Boxplots of decomposition rates by Land Use
plot_data <- model_data %>%
  pivot_longer(cols = c(k_micro, k_total, k_macro),
               names_to = "Decomposition_Type",
               values_to = "k_value")

p1 <- ggplot(plot_data, aes(x = LandUse, y = k_value, fill = LandUse)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  facet_wrap(~ Decomposition_Type, scales = "free_y") +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(title = "Figure 1. Decomposition Rates by Land Use",
       y = "Decomposition Rate (per degree day)",
       x = "Land Use",
       fill = "Land Use") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
print(p1)

# Figure 2: Estimated Marginal Means from Mixed Models
emm_micro <- emmeans(model_micro, ~ LandUse) %>% as.data.frame() %>% mutate(type = "k_micro")
emm_total <- emmeans(model_total, ~ LandUse) %>% as.data.frame() %>% mutate(type = "k_total")
emm_macro <- emmeans(model_macro, ~ LandUse) %>% as.data.frame() %>% mutate(type = "k_macro")

emm_all <- bind_rows(emm_micro, emm_total, emm_macro)

p2 <- ggplot(emm_all, aes(x = LandUse, y = emmean, fill = LandUse)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_wrap(~ type, scales = "free_y") +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(title = "Figure 2. Estimated Marginal Means of Decomposition Rates",
       x = "Land Use",
       y = "Estimated k (per degree day)") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")
print(p2)

# Figure 3: Macroinvertebrate Abundance vs Macroinvertebrate-Mediated Decomposition
p3 <- ggplot(model_data, aes(x = total_abundance, y = k_macro, color = LandUse)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(
    title = "Figure 3. Macroinvertebrate Abundance vs Macroinvertebrate-Mediated Decomposition",
       x = "Total Macroinvertebrate Abundance",
       y = "k_macro (coarse decomposition)",
       color = "Land Use") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        legend.position = "bottom")
print(p3)

# Figure 4: Stream-wise Random Intercepts from Mixed Model
stream_intercepts <- ranef(model_macro)$stream %>%
  as.data.frame() %>%
  rename(deviation = `(Intercept)`) %>%
  rownames_to_column(var = "stream")

p4 <- ggplot(stream_intercepts, aes(x = reorder(stream, deviation), y = deviation)) +
  geom_point(size = 3, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Figure 4. Stream-wise Random Intercepts for k_macro",
       y = "Deviation from Overall Mean",
       x = "Stream")
print(p4)

# Figure 5: Water Temperature vs k_macro
p5 <- ggplot(model_data, aes(x = water_temp, y = k_macro, color = LandUse)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(title = "Figure 5. Water Temperature vs k_macro",
       x = "Water Temperature (°C)",
       y = "k_macro (coarse − fine decomposition)",
       color = "Land Use") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        legend.position = "bottom")
print(p5)

# Figure 6: Partial Effects of Environmental Variables on k_macro
eff_temp <- ggpredict(model_env, terms = "water_temp")
eff_cond <- ggpredict(model_env, terms = "conductivity")
eff_oxy  <- ggpredict(model_env, terms = "oxygen_dissolved")

p6_1 <- plot(eff_temp) + labs(title = "Water Temperature (°C)", y = "Predicted k_macro") + theme_bw()
p6_2 <- plot(eff_cond) + labs(title = "Conductivity (µS/cm)", y = "") + theme_bw()
p6_3 <- plot(eff_oxy) + labs(title = "Dissolved Oxygen (mg/L)", y = "") + theme_bw()

p6 <- p6_1 | p6_2 | p6_3
p6 <- p6 + plot_annotation(title = "Figure 6. Partial Effects of Environmental Predictors on Macroinvertebrate Decomposition")
print(p6)

# ------------------------------------------------------------
# 7. Tables
# ------------------------------------------------------------

get_fixed <- function(model, response_name) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    dplyr::mutate(response = response_name) %>%
    dplyr::select(response, term, estimate, std.error, df, statistic, p.value)
}

# Table 1: Fixed effects summary
fx_micro <- get_fixed(model_micro, "k_micro")
fx_total <- get_fixed(model_total, "k_total")
fx_macro <- get_fixed(model_macro, "k_macro")
fx_macro_abund <- get_fixed(model_macro_abund, "k_macro (+abundance)")

table1 <- bind_rows(fx_micro, fx_total, fx_macro, fx_macro_abund) %>%
  dplyr::mutate(across(c(estimate, std.error, statistic, p.value), ~ round(.x, 3)))
print(table1)

# Fit linear model for each LandUse group (example: k_macro ~ total_abundance)
library(dplyr)
library(ggplot2)

# Plot measured means and model fit, faceted by LandUse
p_ribbon <- ggplot(model_data, aes(x = total_abundance, y = k_macro, color = LandUse)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = LandUse), alpha = 0.25, linetype = "solid") +
  facet_wrap(~LandUse, scales = "free") +
  labs(
    title = "Macroinvertebrate-Mediated Decomposition vs. Invertebrate Abundance",
    x = "Total Macroinvertebrate Abundance",
    y = "k_macro (coarse − fine decomposition)"
  ) +
  theme_minimal()
print(p_ribbon)

# Facet by stream
ggplot(model_data, aes(x = total_abundance, y = k_macro, color = LandUse)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = LandUse), alpha = 0.25) +
  facet_wrap(~stream, scales = "free") +
  labs(
    title = "k_macro vs. Abundance, by Stream",
    x = "Total Macroinvertebrate Abundance",
    y = "k_macro"
  ) +
  theme_minimal()


summary_table <- model_data %>%
  group_by(LandUse, stream) %>%
  summarise(
    n = n(),
    mean_k_macro = mean(k_macro, na.rm = TRUE),
    sd_k_macro = sd(k_macro, na.rm = TRUE)
  )
print(summary_table)



# Table 3: Pairwise contrasts
cont_micro <- emmeans(model_micro, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_micro")

cont_total <- emmeans(model_total, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_total")

cont_macro <- emmeans(model_macro, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_macro")

cont_macro_abund <- emmeans(model_macro_abund, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_macro (+abundance)")

table3 <- bind_rows(cont_micro, cont_total, cont_macro, cont_macro_abund) %>%
  dplyr::select(response, contrast, estimate, SE, df, t.ratio, p.value) %>%
  dplyr::mutate(across(c(estimate, SE, t.ratio, p.value), ~ round(.x, 3)))
print(table3)


# Calculate macroinvertebrate abundance and percent composition per stream
macro_abund_table <- inverts_clean %>%
  group_by(stream) %>%
  summarise(
    Gammarids_n    = sum(Gammarids, na.rm = TRUE),
    Trichoptera_n  = sum(Trichoptera, na.rm = TRUE),
    Chironomids_n  = sum(Chironomids, na.rm = TRUE),
    Others_n       = sum(Others, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    Total_n = Gammarids_n + Trichoptera_n + Chironomids_n + Others_n,
    Gammarids_pct   = paste0(round(Gammarids_n / Total_n * 100, 1), "%"),
    Trichoptera_pct = paste0(round(Trichoptera_n / Total_n * 100, 1), "%"),
    Chironomids_pct = paste0(round(Chironomids_n / Total_n * 100, 1), "%"),
    Others_pct      = paste0(round(Others_n / Total_n * 100, 1), "%")
  ) %>%
  dplyr::select(
    stream,
    Gammarids_n,   Gammarids_pct,
    Trichoptera_n, Trichoptera_pct,
    Chironomids_n, Chironomids_pct,
    Others_n,      Others_pct,
    Total_n
  )

# View the table
print(macro_abund_table)

# Merge LandUse info into inverts_clean for grouping
inverts_LandUse <- inverts_clean %>%
  left_join(model_data %>% select(stream, LandUse) %>% distinct(), by = "stream")

# Calculate totals and percent for each taxon per LandUse group
macro_abund_landuse <- inverts_LandUse %>%
  dplyr::filter(!is.na(LandUse)) %>%
  group_by(LandUse) %>%
  summarise(
    Gammarids_n    = sum(Gammarids, na.rm = TRUE),
    Trichoptera_n  = sum(Trichoptera, na.rm = TRUE),
    Chironomids_n  = sum(Chironomids, na.rm = TRUE),
    Others_n       = sum(Others, na.rm = TRUE)
  ) %>%
  mutate(
    Total_n = Gammarids_n + Trichoptera_n + Chironomids_n + Others_n,
    Gammarids_pct   = paste0(round(Gammarids_n / Total_n * 100, 1), "%"),
    Trichoptera_pct = paste0(round(Trichoptera_n / Total_n * 100, 1), "%"),
    Chironomids_pct = paste0(round(Chironomids_n / Total_n * 100, 1), "%"),
    Others_pct      = paste0(round(Others_n / Total_n * 100, 1), "%")
  ) %>%
  dplyr::select(
    LandUse,
    Gammarids_n,   Gammarids_pct,
    Trichoptera_n, Trichoptera_pct,
    Chironomids_n, Chironomids_pct,
    Others_n,      Others_pct,
    Total_n
  )
# Allow LandUse effect to vary by stream
model_micro_slope <- lmer(k_micro ~ LandUse + (LandUse | stream), data = model_data)
summary(model_micro_slope)


# Show the result
print(macro_abund_landuse)

# For RMarkdown, print as a table:
knitr::kable(macro_abund_landuse, caption = "Macroinvertebrate Abundance and Percent by Land Use")


