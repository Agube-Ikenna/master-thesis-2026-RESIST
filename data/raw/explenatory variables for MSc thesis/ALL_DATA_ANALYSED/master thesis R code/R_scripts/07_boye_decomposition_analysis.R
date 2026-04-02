# ==================================================================
# Leaf Decomposition in the Boye Catchment – Full R Workflow
# Author: Agube Ikenna Stephen
# Supervisor: Dr. Verena C. Schreiner
# Department of Ecotoxicology, University of Duisburg-Essen
# Description: Mixed model analysis – restoration status as predictor
# ==================================================================

# NOTE: The Boye analysis differs fundamentally from the Kinzig analysis.
# Kinzig compares LAND USE (rural vs urban).
# Boye compares RESTORATION STATUS:
#   - reference: BOYkl (Kleine Boye, unrestored reference)
#   - recently_restored: BOYohB224, BOYohKi, BOYuhHa (restored 2021, ~1 yr recovery)
#   - long_restored: BOYohBr, BOYohSp, BOYuhSp (restored 2002, ~20 yr recovery)
#
# Research question: Does leaf litter decomposition recover with stream restoration,
# and how does recovery time affect microbial vs macroinvertebrate-mediated breakdown?

# ------------------------------------------------------------
# 1. Setup: Load Libraries and Set Working Directory
# ------------------------------------------------------------

# ── Working directory (portable – no manual editing needed) ───────────────────
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  script_dir <- getSrcDirectory(function() {})
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

# These files contain data for ALL catchments (Boye + Kinzig).
# We filter to Boye below.
processed    <- read_excel("2025_field_processed_samples.xlsx")
initial_mass <- read_excel("2025_field_initial_mass.xlsx")
transport    <- read_excel("2025_field_transport_controls.xlsx")
inverts      <- read_excel("2025_field_invertebrates_coarse_leaf_bags.xlsx")

# Boye-specific physicochemistry (different file from Kinzig, different column names)
phch_raw <- read_excel("boye_processed_samples_mar2025.xlsx") %>% clean_names()

# Check column names of phch file to confirm the join key
# names(phch_raw)

# ------------------------------------------------------------
# 3. Define Restoration Status for Each Boye Site
# ------------------------------------------------------------

# Source: SFB_BasicData_SamplingSites_Boye_&_Kinzig file (see data/processed/)
# BOYkl          = Kleine Boye (reference, no restoration)
# BOYohB224      = restored 2021 (recently restored, 1 yr recovery)
# BOYohBr        = restored 2002 (long-restored, ~20 yr recovery)
# BOYohKi        = restored 2021 (recently restored, 1 yr recovery)
# BOYohSp        = restored 2002 (long-restored, ~20 yr recovery)
# BOYuhHa        = restored 2021 (recently restored, 1 yr recovery)
# BOYuhSp        = restored 2002 (long-restored, ~20 yr recovery)
#
# NOTE: stream names in processed data are mixed case (e.g., "BOYkl", "BOYohB224").
# Verify with: unique(processed$stream[processed$catchment == "Boye"])

restoration_status <- tibble(
  stream = c("BOYkl", "BOYohB224", "BOYohBr", "BOYohKi", "BOYohSp", "BOYuhHa", "BOYuhSp"),
  RestStatus = c("reference",
                 "recently_restored",
                 "long_restored",
                 "recently_restored",
                 "long_restored",
                 "recently_restored",
                 "long_restored")
)

# Make RestStatus an ordered factor so models reflect the gradient of recovery
restoration_status <- restoration_status %>%
  mutate(RestStatus = factor(RestStatus,
                             levels = c("reference", "recently_restored", "long_restored"),
                             ordered = FALSE))

# ------------------------------------------------------------
# 4. Calculate Decomposition Rates (k values)
# ------------------------------------------------------------

# Handling loss correction from transport controls
transport <- transport %>% mutate(AFDM_0 = mass_60_g - mass_500_g)
mean_AFDM_0 <- mean(transport$AFDM_0, na.rm = TRUE)

# Filter processed data to Boye only and compute k per degree day
boye_processed <- processed %>%
  dplyr::filter(catchment == "Boye") %>%
  dplyr::mutate(
    OMD_60_mg  = as.numeric(OMD_60_mg),
    OMD_500_mg = as.numeric(OMD_500_mg),
    AFDM_final = OMD_60_mg - OMD_500_mg,
    AFDM_0     = mean_AFDM_0,
    # NOTE: Update degree-days below if Boye deployment duration or temperature differs from Kinzig
    # Kinzig used: 14 days × 12.5 °C = 175 degree-days
    # Check boye_temperature_logger_march2025.csv for actual mean temperature at Boye
    deg_days   = 14 * 12.5,
    OMD_dday   = ((AFDM_0 - AFDM_final) / AFDM_0) * 100 / deg_days
  )

# Verify which stream names appear in the Boye data:
message("Boye streams found in processed data: ")
print(unique(boye_processed$stream))

# Merge restoration status
boye_labeled <- boye_processed %>%
  left_join(restoration_status, by = "stream") %>%
  filter(!is.na(RestStatus), !is.na(OMD_dday))

# Reshape to wide format: fine = microbial, coarse = total, difference = macroinvertebrate
model_data <- boye_labeled %>%
  dplyr::select(stream, replicate, type, RestStatus, OMD_dday) %>%
  pivot_wider(names_from = type, values_from = OMD_dday) %>%
  rename(k_micro = fine, k_total = coarse) %>%
  dplyr::mutate(k_macro = k_total - k_micro)

# ------------------------------------------------------------
# 5. Join Invertebrate Abundance Data (Coarse Leaf Bags)
# ------------------------------------------------------------

inverts_clean <- inverts %>%
  dplyr::filter(catchment == "Boye", type == "coarse") %>%
  dplyr::mutate(across(c(Gammarids, Trichoptera, Chironomids, Others), as.numeric)) %>%
  dplyr::mutate(total_abundance = rowSums(across(c(Gammarids, Trichoptera, Chironomids, Others)),
                                          na.rm = TRUE)) %>%
  dplyr::select(stream, replicate, total_abundance, Gammarids, Trichoptera, Chironomids, Others)

model_data <- model_data %>%
  dplyr::mutate(across(c(stream, replicate), str_trim)) %>%
  left_join(inverts_clean, by = c("stream", "replicate"))

# ------------------------------------------------------------
# 6. Join Physicochemical Data
# ------------------------------------------------------------

# NOTE: The Boye phch file uses "site" as the join key, not "stream".
# Column names after clean_names(): site, p_h_dimensionless, conductivity_u_s_cm,
# dissolved_oxygen_mg_l, oxygen_saturation_pct, water_temp_deg_c, water_level_cm, etc.
# Check with: names(phch_raw) before running this block.

phch_clean <- phch_raw %>%
  rename(stream = site) %>%                        # rename join key to match model_data
  dplyr::select(
    stream,
    ph              = p_h_dimensionless,
    water_temp      = water_temp_deg_c,
    oxygen_dissolved = dissolved_oxygen_mg_l,
    oxygen_saturation = oxygen_saturation_pct,
    conductivity    = conductivity_u_s_cm,
    water_level     = water_level_cm
  ) %>%
  dplyr::mutate(stream = str_trim(stream))

model_data <- model_data %>%
  left_join(phch_clean, by = "stream")

# Check the merged structure
glimpse(model_data)
message("Rows in final model_data: ", nrow(model_data))
message("Sites covered: ", paste(unique(model_data$stream), collapse = ", "))

# ------------------------------------------------------------
# 7. Fit Linear Mixed-Effects Models
# ------------------------------------------------------------

# HYPOTHESIS 1: Microbial decomposition differs by restoration status
model_micro <- lmer(k_micro ~ RestStatus + (1 | stream), data = model_data)
summary(model_micro)

# HYPOTHESIS 2: Total decomposition differs by restoration status
model_total <- lmer(k_total ~ RestStatus + (1 | stream), data = model_data)
summary(model_total)

# HYPOTHESIS 3: Macroinvertebrate-mediated decomposition differs by restoration status
model_macro <- lmer(k_macro ~ RestStatus + (1 | stream), data = model_data)
summary(model_macro)

# Extended: add invertebrate abundance as covariate
model_macro_abund <- lmer(k_macro ~ RestStatus + total_abundance + (1 | stream), data = model_data)
summary(model_macro_abund)

# Environmental covariates model
model_env <- lmer(k_macro ~ RestStatus + water_temp + conductivity + oxygen_dissolved + (1 | stream),
                  data = model_data)
summary(model_env)

# ------------------------------------------------------------
# 8. Figures
# ------------------------------------------------------------

# Figure 1: Boxplots of decomposition rates by Restoration Status
plot_data <- model_data %>%
  pivot_longer(cols = c(k_micro, k_total, k_macro),
               names_to  = "Decomposition_Type",
               values_to = "k_value")

p1 <- ggplot(plot_data, aes(x = RestStatus, y = k_value, fill = RestStatus)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  facet_wrap(~ Decomposition_Type, scales = "free_y") +
  scale_fill_manual(values = c("reference"          = "gray50",
                                "recently_restored"  = "steelblue",
                                "long_restored"      = "forestgreen")) +
  labs(title = "Figure 1. Decomposition Rates by Restoration Status (Boye)",
       y     = "Decomposition Rate (per degree day)",
       x     = "Restoration Status",
       fill  = "Restoration Status") +
  theme_bw() +
  theme(plot.title   = element_text(size = 14, face = "bold"),
        strip.text   = element_text(face = "bold"),
        axis.text.x  = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
print(p1)

# Figure 2: Estimated Marginal Means
emm_micro <- emmeans(model_micro, ~ RestStatus) %>% as.data.frame() %>% mutate(type = "k_micro")
emm_total <- emmeans(model_total, ~ RestStatus) %>% as.data.frame() %>% mutate(type = "k_total")
emm_macro <- emmeans(model_macro, ~ RestStatus) %>% as.data.frame() %>% mutate(type = "k_macro")

emm_all <- bind_rows(emm_micro, emm_total, emm_macro)

p2 <- ggplot(emm_all, aes(x = RestStatus, y = emmean, fill = RestStatus)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_wrap(~ type, scales = "free_y") +
  scale_fill_manual(values = c("reference"          = "gray50",
                                "recently_restored"  = "steelblue",
                                "long_restored"      = "forestgreen")) +
  labs(title = "Figure 2. Estimated Marginal Means – Boye Restoration Sites",
       x     = "Restoration Status",
       y     = "Estimated k (per degree day)") +
  theme_bw() +
  theme(plot.title  = element_text(size = 14, face = "bold"),
        strip.text  = element_text(face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none")
print(p2)

# Figure 3: Invertebrate Abundance vs Macroinvertebrate-Mediated Decomposition
p3 <- ggplot(model_data, aes(x = total_abundance, y = k_macro, color = RestStatus)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = RestStatus), alpha = 0.2) +
  facet_wrap(~ RestStatus, scales = "free") +
  scale_color_manual(values = c("reference"         = "gray50",
                                "recently_restored" = "steelblue",
                                "long_restored"     = "forestgreen")) +
  scale_fill_manual(values  = c("reference"         = "gray50",
                                "recently_restored" = "steelblue",
                                "long_restored"     = "forestgreen")) +
  labs(title = "Figure 3. Invertebrate Abundance vs k_macro (Boye)",
       x     = "Total Macroinvertebrate Abundance",
       y     = "k_macro (coarse – fine decomposition)") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")
print(p3)

# Figure 4: Stream-wise Random Intercepts
stream_intercepts <- ranef(model_macro)$stream %>%
  as.data.frame() %>%
  rename(deviation = `(Intercept)`) %>%
  rownames_to_column(var = "stream") %>%
  left_join(restoration_status, by = "stream")

p4 <- ggplot(stream_intercepts, aes(x = reorder(stream, deviation), y = deviation,
                                     color = RestStatus)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  coord_flip() +
  scale_color_manual(values = c("reference"         = "gray50",
                                "recently_restored" = "steelblue",
                                "long_restored"     = "forestgreen")) +
  labs(title = "Figure 4. Stream-wise Random Intercepts for k_macro (Boye)",
       y     = "Deviation from Overall Mean",
       x     = "Stream") +
  theme_minimal()
print(p4)

# Figure 5: Partial Effects of Environmental Predictors
eff_temp <- ggpredict(model_env, terms = "water_temp")
eff_cond <- ggpredict(model_env, terms = "conductivity")
eff_oxy  <- ggpredict(model_env, terms = "oxygen_dissolved")

p5_1 <- plot(eff_temp) + labs(title = "Water Temperature (°C)", y = "Predicted k_macro") + theme_bw()
p5_2 <- plot(eff_cond) + labs(title = "Conductivity (µS/cm)", y = "") + theme_bw()
p5_3 <- plot(eff_oxy)  + labs(title = "Dissolved Oxygen (mg/L)", y = "") + theme_bw()

p5 <- (p5_1 | p5_2 | p5_3) +
  plot_annotation(title = "Figure 5. Partial Effects of Environmental Predictors on k_macro (Boye)")
print(p5)

# ------------------------------------------------------------
# 9. Tables
# ------------------------------------------------------------

get_fixed <- function(model, response_name) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    dplyr::mutate(response = response_name) %>%
    dplyr::select(response, term, estimate, std.error, df, statistic, p.value)
}

fx_micro       <- get_fixed(model_micro,       "k_micro")
fx_total       <- get_fixed(model_total,       "k_total")
fx_macro       <- get_fixed(model_macro,       "k_macro")
fx_macro_abund <- get_fixed(model_macro_abund, "k_macro (+abundance)")

table1 <- bind_rows(fx_micro, fx_total, fx_macro, fx_macro_abund) %>%
  dplyr::mutate(across(c(estimate, std.error, statistic, p.value), ~ round(.x, 3)))
print(table1)

# Pairwise contrasts between restoration status groups
cont_micro <- emmeans(model_micro, pairwise ~ RestStatus)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_micro")
cont_total <- emmeans(model_total, pairwise ~ RestStatus)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_total")
cont_macro <- emmeans(model_macro, pairwise ~ RestStatus)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_macro")

table2 <- bind_rows(cont_micro, cont_total, cont_macro) %>%
  dplyr::select(response, contrast, estimate, SE, df, t.ratio, p.value) %>%
  dplyr::mutate(across(c(estimate, SE, t.ratio, p.value), ~ round(.x, 3)))
print(table2)

# Summary table: mean k values per site and restoration status
summary_table <- model_data %>%
  group_by(RestStatus, stream) %>%
  summarise(
    n           = n(),
    mean_k_micro = mean(k_micro, na.rm = TRUE),
    mean_k_total = mean(k_total, na.rm = TRUE),
    mean_k_macro = mean(k_macro, na.rm = TRUE),
    sd_k_macro   = sd(k_macro, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_table)
