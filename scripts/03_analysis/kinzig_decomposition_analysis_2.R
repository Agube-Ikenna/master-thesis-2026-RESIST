# ==========================================================================
# The Effect of Multiple Stressors on Microbial and Total Leaf Decomposition
# in Urban and Rural Streams of the Kinzig Catchment
#
# Author: Ikenna Stephen Agube
# Supervisor: Dr. Verena C. Schreiner
# Department of Ecotoxicology, University of Duisburg-Essen
# Date: 23.07.2025
# Contact: ikenna.agube@uni-due.de; verena.schreiner@rptu.de
# ==========================================================================

# Description:
#   This script implements a clean and reproducible workflow to analyse how
#   leaf decomposition rates differ between rural and urban streams in the
#   Kinzig catchment. It uses properly corrected initial masses, accounts for
#   handling losses via transport controls, corrects for disks removed for
#   other analyses, and calculates the exponential decomposition rate (k)
#   normalised by degree days. Mixed‑effects models evaluate differences in
#   microbial, total and macroinvertebrate‑mediated decomposition between
#   land uses, while summary statistics and simple plots communicate the
#   results. Each step is documented with comments so that future users can
#   understand and reproduce the analysis.
# ============================================================================

## --------------------------------------------------------------------------
## 1) Load required packages and set general options
## --------------------------------------------------------------------------

library(tidyverse)      # core data manipulation and plotting
library(readxl)         # read Excel files (not used here but loaded for completeness)
library(lme4)           # fit mixed‑effects models
library(lmerTest)       # obtain p‑values for lme4 models
library(performance)    # check model assumptions
library(janitor)        # clean column names
library(ggeffects)      # generate predicted values and partial effects
library(patchwork)      # arrange ggplots
library(broom.mixed)    # convert mixed models to tidy data frames
library(stringr)        # string manipulation
library(glue)           # formatted output for reporting

# Improve readability of numeric output (avoid scientific notation)
options(scipen = 999, digits = 3)

## --------------------------------------------------------------------------
## 2) Define working directory and import data
## --------------------------------------------------------------------------

# Please, This is the only line that needs to change when moving to a different machine.

setwd("C:/Users/agube/OneDrive/Documentos/r_programming_2025")
message(glue("Working directory set to: {getwd()}"))

# When reading files exported from Excel, specify NA strings so that
# missing values are treated consistently. Using CSV instead of XLSX
# avoids the subtle floating point changes observed in Excel imports.

processed    <- read.csv("2025_field_processed_samples.csv", na.strings = c("NR", "NA"))
initial_mass <- read.csv("2025_field_initial_mass.csv",    na.strings = c("NR", "NA"))
transport    <- read.csv("2025_field_transport_controls.csv", na.strings = c("NR", "NA"))
site_info    <- read.csv("2025_field_site_information_Kinzig.csv", na.strings = c("NR", "NA"))
inverts      <- read.csv("2025_field_invertebrates_coarse_leaf_bags.csv", na.strings = c("NR", "NA"))
phch_raw     <- read.csv("202505_phch_data_Kinzig.csv", na.strings = c("NR", "NA"))

## --------------------------------------------------------------------------
## 3) Clean basic variables: dates, numeric types and land use
## --------------------------------------------------------------------------

# Convert deployment and sampling dates to Date objects and compute
# deployment duration in days. Converting early allows us to use the
# difference when calculating degree days later.
processed$date_deployment <- as.Date(processed$date_deployment, format = "%d/%m/%Y")
processed$date_sampled    <- as.Date(processed$date_sampled,    format = "%d/%m/%Y")
processed$deployment_duration_days <- as.numeric(difftime(processed$date_sampled,
                                                          processed$date_deployment,
                                                          units = "days"))

# Filter to Kinzig catchment only. If you wish to run analyses on another
# catchment later, adapt this line accordingly.
processed_kinzig <- processed %>% filter(catchment == "Kinzig")

# Ensure that numeric columns are indeed numeric. Files with missing values
# often import numeric fields as characters; this loop coerces them safely.
numeric_cols <- c("OMD_60_mg", "OMD_500_mg", "no_disks_OMD", "no_disks_all",
                  "rest_60_g", "rest_500_g", "initial_mass_g")
for (col in numeric_cols) {
  if (col %in% names(processed_kinzig)) {
    processed_kinzig[[col]] <- as.numeric(processed_kinzig[[col]])
  }
}

# Translate the land use flags from site_info into a factor with levels
# "rural" and "urban". Any stream lacking a land use flag becomes NA and is
# dropped later. We keep only the stream and LandUse columns from site_info.
site_info <- site_info %>%
  mutate(LandUse = case_when(
    tolower(str_trim(rural)) == "x" ~ "rural",
    tolower(str_trim(urban)) == "x" ~ "urban",
    TRUE ~ NA_character_
  )) %>%
  select(stream, LandUse)

processed_kinzig <- processed_kinzig %>%
  left_join(site_info, by = "stream")

## --------------------------------------------------------------------------
## 4) Handling‑loss correction and AFDM calculations
## --------------------------------------------------------------------------

# To correct the initial mass for handling losses, we calculate the
# fraction of ash‑free dry mass (AFDM) lost between 60 °C and 500 °C in
# the transport controls. Because handling losses differ by mesh type
# (fine versus coarse), we compute the mean AFDM fraction separately for
# each type. The transport file uses lowercase catchment and type labels
# for consistency.

transport <- transport %>%
  mutate(catchment = tolower(str_trim(catchment)),
         type      = tolower(str_trim(type))) %>%
  mutate(
    perc_60   = mass_60_g  / initial_mass_g,
    perc_500  = mass_500_g / initial_mass_g,
    afdm_frac = perc_60 - perc_500
  ) %>%
  filter(is.finite(afdm_frac), afdm_frac > 0, afdm_frac < 1)

# Compute the mean AFDM fraction for Kinzig by mesh type. If the
# transport controls contain other catchments, they are ignored. We
# require at least one control per type; otherwise, the analysis would
# stop.
corr_afdm_kinzig <- transport %>%
  filter(catchment == "kinzig") %>%
  group_by(type) %>%
  summarise(afdm_frac = mean(afdm_frac, na.rm = TRUE), .groups = "drop")
stopifnot(nrow(corr_afdm_kinzig) >= 1)

# Bring in initial masses if missing. Some rows in processed_kinzig may
# lack initial_mass_g because the values reside in a separate file. We
# match by bag_ID and type.
if (!"initial_mass_g" %in% names(processed_kinzig)) {
  initial_mass <- initial_mass %>% mutate(type = tolower(str_trim(type)))
  processed_kinzig <- processed_kinzig %>%
    left_join(initial_mass %>% select(bag_ID, type, initial_mass_g), by = c("bag_ID", "type"))
}

# Standardise type labels to lower case for joins
processed_kinzig <- processed_kinzig %>% mutate(type = tolower(str_trim(type)))

# Compute AFDM at the start using the mean handling loss fraction.
processed_kinzig <- processed_kinzig %>%
  left_join(corr_afdm_kinzig, by = "type") %>%
  mutate(AFDM_start = initial_mass_g * afdm_frac) %>%
  relocate(AFDM_start, .after = initial_mass_g)

# Correct the final AFDM for disks removed to measure other endpoints. We
# calculate the mean AFDM per disk for each bag; if unavailable we fall
# back to the stream mean and then to the global mean. The conversion
# from mg (OMD measurements) to g occurs at the end.
processed_kinzig <- processed_kinzig %>%
  mutate(no_disks_OMD = replace_na(no_disks_OMD, 0),
         no_disks_all = replace_na(no_disks_all, 0)) %>%
  mutate(
    mean_afdm_per_disk = ifelse(is.finite(OMD_60_mg) & is.finite(OMD_500_mg) & no_disks_OMD > 0,
                                (OMD_60_mg - OMD_500_mg) / no_disks_OMD,
                                NA_real_)
  )

# Stream‑level mean AFDM per disk
stream_means <- processed_kinzig %>%
  group_by(stream) %>%
  summarise(stream_mean_disk = mean(mean_afdm_per_disk, na.rm = TRUE), .groups = "drop")

# Global mean AFDM per disk (in mg)
global_mean <- processed_kinzig %>%
  filter(is.finite(mean_afdm_per_disk)) %>%
  summarise(val = mean(mean_afdm_per_disk, na.rm = TRUE)) %>%
  pull(val)
if (!is.finite(global_mean)) global_mean <- 0

processed_kinzig <- processed_kinzig %>%
  left_join(stream_means, by = "stream") %>%
  mutate(
    mean_afdm_per_disk = coalesce(mean_afdm_per_disk, stream_mean_disk, global_mean),
    afdm_all_disks     = (no_disks_all * mean_afdm_per_disk) / 1000,  # mg to g
    rest_AFDM          = rest_60_g - rest_500_g,
    AFDM_end           = rest_AFDM + afdm_all_disks
  ) %>%
  select(-stream_mean_disk)

# Remove any rows where core values are missing; models cannot handle NAs.
n0 <- nrow(processed_kinzig)
processed_kinzig <- processed_kinzig %>%
  filter(is.finite(rest_60_g), is.finite(rest_500_g), is.finite(AFDM_start), is.finite(AFDM_end))
message(glue("Dropped {n0 - nrow(processed_kinzig)} rows due to missing masses."))

## --------------------------------------------------------------------------
## 5) Degree days and exponential decomposition rate k
## --------------------------------------------------------------------------

# Summarise physico‑chemical measurements by stream; we keep the mean
# water temperature for computing degree days. When multiple
# observations exist per stream, we take their average.
phch_summ <- phch_raw %>%
  clean_names() %>%
  group_by(stream) %>%
  summarise(water_temp = mean(water_temp, na.rm = TRUE), .groups = "drop")

processed_kinzig <- processed_kinzig %>%
  left_join(phch_summ, by = "stream") %>%
  mutate(
    deg_days = deployment_duration_days * water_temp,
    k        = -log(AFDM_end / AFDM_start) / deg_days
  ) %>%
  filter(is.finite(k), k > 0, AFDM_end < AFDM_start, deg_days > 0) %>%
  mutate(
    # round key variables for reporting (not needed for computation)
    AFDM_start     = round(AFDM_start, 4),
    AFDM_end       = round(AFDM_end, 4),
    rest_AFDM      = round(rest_AFDM, 4),
    afdm_all_disks = round(afdm_all_disks, 4),
    k              = round(k, 5)
  )

# Drop rows without a land use assignment
processed_kinzig <- processed_kinzig %>% filter(!is.na(LandUse))

# Set factor levels explicitly: rural as the reference level
processed_kinzig$LandUse <- factor(processed_kinzig$LandUse, levels = c("rural", "urban"))

## --------------------------------------------------------------------------
## 6) Prepare datasets for hypotheses H1–H3
## --------------------------------------------------------------------------

# Separate microbial (fine mesh) and total (coarse mesh) data. We keep
# only rows with finite k > 0. Later we pair fine and coarse data to
# compute macroinvertebrate‑mediated decomposition.
microbial_data <- processed_kinzig %>% filter(type == "fine"   & is.finite(k) & k > 0)
total_data     <- processed_kinzig %>% filter(type == "coarse" & is.finite(k) & k > 0)

# Pair coarse and fine bags by stream, replicate and land use to compute
# macroinvertebrate contribution. We take care to subtract within
# matching bags only. Cases where the difference is non‑positive are
# removed, as negative contributions are biologically implausible or
# indicate measurement error.
microbial_matched <- microbial_data %>%
  select(stream, replicate, LandUse, k_microbial = k) %>%
  left_join(total_data %>% select(stream, replicate, LandUse, k_total = k),
            by = c("stream", "replicate", "LandUse")) %>%
  mutate(k_macro = k_total - k_microbial) %>%
  filter(is.finite(k_macro), k_macro > 0)

message(glue("N fine bags: {nrow(microbial_data)}, N coarse bags: {nrow(total_data)}, N paired macro bags: {nrow(microbial_matched)}"))

## --------------------------------------------------------------------------
## 7) Fit mixed‑effects models for H1–H3
## --------------------------------------------------------------------------

# We model k (or k_macro) as a function of LandUse with a random
# intercept for stream. lmerTest supplies p‑values for fixed effects.

model_h1_macro <- lmer(k_macro ~ LandUse + (1 | stream), data = microbial_matched)
model_h2_micro <- lmer(k ~ LandUse + (1 | stream), data = microbial_data)
model_h3_total <- lmer(k ~ LandUse + (1 | stream), data = total_data)

# Perform diagnostic checks (residuals, QQ‑plot, etc.). These return
# ggplot objects. Inspect them interactively to evaluate whether model
# assumptions are violated. We do not include them in the final report
# automatically, but encourage you to check them.
check_h1 <- performance::check_model(model_h1_macro)
check_h2 <- performance::check_model(model_h2_micro)
check_h3 <- performance::check_model(model_h3_total)

# Extract tidy summaries for the fixed effects. In these models
# LandUseurban gives the difference (urban – rural). A negative
# estimate implies decomposition is lower in urban streams.
tidy_h1_macro <- tidy(model_h1_macro, effects = "fixed")
tidy_h2_micro <- tidy(model_h2_micro, effects = "fixed")
tidy_h3_total <- tidy(model_h3_total, effects = "fixed")

extract_effect <- function(tdf) {
  tibble(
    est = tdf %>% filter(term == "LandUseurban") %>% pull(estimate),
    se  = tdf %>% filter(term == "LandUseurban") %>% pull(std.error),
    p   = tdf %>% filter(term == "LandUseurban") %>% pull(p.value)
  )
}

eff_h1 <- extract_effect(tidy_h1_macro)
eff_h2 <- extract_effect(tidy_h2_micro)
eff_h3 <- extract_effect(tidy_h3_total)
message(glue("H1 (k_macro)  urban − rural = {round(eff_h1$est,5)} ± {round(eff_h1$se,5)}; p = {signif(eff_h1$p,3)}"))
message(glue("H2 (k_micro)  urban − rural = {round(eff_h2$est,5)} ± {round(eff_h2$se,5)}; p = {signif(eff_h2$p,3)}"))
message(glue("H3 (k_total)  urban − rural = {round(eff_h3$est,5)} ± {round(eff_h3$se,5)}; p = {signif(eff_h3$p,3)}"))

## --------------------------------------------------------------------------
## 8) Optional: Model with abundance as covariate (context only)
## --------------------------------------------------------------------------

# Although not part of the three hypotheses, it can be biologically
# interesting to see whether the total abundance of invertebrates
# correlates with k. We merge total abundance into the main dataset and
# fit a model including abundance as an additional fixed effect. We do
# not interpret temperature again because degree days already incorporate
# temperature. Only inspect the sign of the abundance coefficient.
inverts_clean <- inverts %>%
  mutate(across(c(Gammarids, Trichoptera, Chironomids, Others), as.numeric)) %>%
  mutate(total_abundance = rowSums(across(c(Gammarids, Trichoptera, Chironomids, Others)), na.rm = TRUE)) %>%
  select(stream, replicate, total_abundance)

processed_kinzig_abund <- processed_kinzig %>%
  left_join(inverts_clean, by = c("stream", "replicate"))

model_k_abund <- lmer(k ~ LandUse + total_abundance + (1 | stream), data = processed_kinzig_abund)
tidy_k_abund <- tidy(model_k_abund, effects = "fixed")

## --------------------------------------------------------------------------
## 9) Generate figures to visualise the results
## --------------------------------------------------------------------------

# Boxplots summarising microbial, total and macro rates by land use. A
# panel of three plots helps the reader compare the distributions.
p_microbial <- ggplot(microbial_data, aes(x = LandUse, y = k, fill = LandUse)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(title = "Microbial decomposition rate", x = "Land use", y = "k (per degree day)") +
  theme_bw() + theme(legend.position = "none")

p_total <- ggplot(total_data, aes(x = LandUse, y = k, fill = LandUse)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(title = "Total decomposition rate", x = "Land use", y = "k (per degree day)") +
  theme_bw() + theme(legend.position = "none")

p_macro <- ggplot(microbial_matched, aes(x = LandUse, y = k_macro, fill = LandUse)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(title = "Macroinvertebrate contribution", x = "Land use", y = "k_macro (per degree day)") +
  theme_bw() + theme(legend.position = "none")

final_panel <- p_microbial | p_total | p_macro

print(final_panel + plot_annotation(
  caption = "Microbial: fine mesh; Total: coarse mesh; Macro: k_total − k_micro"
))

## --------------------------------------------------------------------------
## 10) Summarise results in a compact table
## --------------------------------------------------------------------------

# Means and standard errors by land use for each endpoint
summary_all <- bind_rows(
  microbial_data %>% group_by(LandUse) %>%
    summarise(type = "Microbial", n = n(),
              mean_k = mean(k, na.rm = TRUE),
              se_k   = sd(k,  na.rm = TRUE)/sqrt(n()), .groups = "drop"),
  total_data %>% group_by(LandUse) %>%
    summarise(type = "Total", n = n(),
              mean_k = mean(k, na.rm = TRUE),
              se_k   = sd(k,  na.rm = TRUE)/sqrt(n()), .groups = "drop"),
  microbial_matched %>% group_by(LandUse) %>%
    summarise(type = "Macroinvertebrate", n = n(),
              mean_k = mean(k_macro, na.rm = TRUE),
              se_k   = sd(k_macro, na.rm = TRUE)/sqrt(n()), .groups = "drop")
)

# Combine fixed effect summaries into one table for reporting. The row
# corresponding to "LandUseurban" shows the urban–rural difference.
table_models <- bind_rows(
  tidy_h1_macro %>% mutate(Endpoint = "Macroinvertebrate"),
  tidy_h2_micro %>% mutate(Endpoint = "Microbial"),
  tidy_h3_total %>% mutate(Endpoint = "Total")
)
# export_figures_tables.R
# Save combined panel of decomposition rates
png("decomposition_rates_panel.png", width = 9, height = 3, units = "in", res = 300)
print(final_panel + plot_annotation(
  caption = "Microbial: fine mesh; Total: coarse mesh; Macro: k_total − k_micro"
))
dev.off()

# Optionally save individual boxplots
ggsave("microbial_boxplot.png", p_microbial, width = 3, height = 4, dpi = 300)
ggsave("total_boxplot.png",     p_total,     width = 3, height = 4, dpi = 300)
ggsave("macro_boxplot.png",     p_macro,     width = 3, height = 4, dpi = 300)

# Save summary statistics by land use and endpoint
write.csv(summary_all, "decomposition_summary_table.csv", row.names = FALSE)

# Save fixed-effects table from mixed models
write.csv(table_models, "mixed_model_fixed_effects_table.csv", row.names = FALSE)

# Save microbial, total and macro plots separately
ggsave("microbial_boxplot.png", p_microbial, width = 3, height = 4, dpi = 300)
ggsave("total_boxplot.png",     p_total,     width = 3, height = 4, dpi = 300)
ggsave("macro_boxplot.png",     p_macro,     width = 3, height = 4, dpi = 300)

# Print the summary tables. In an RMarkdown context, use knitr::kable()
# to produce nicely formatted tables; here we simply display the
# data frames.
print(summary_all)
print(table_models)

# End of script
