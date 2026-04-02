##############################################################################
# LEAF DECOMPOSITION STUDY — COMPLETE R ANALYSIS SCRIPT
# Based on Schreiner et al. (2023), Environmental Toxicology and Chemistry
# Study streams: Boye (urban/industrial) | Kinzig (agricultural/rural)
# Data collection: March 2025 | Field deployment ~15 days
#
# Author: [Your Name] | Date: March 2025
# Contact: [your.email@institution.de]
##############################################################################

# ============================================================================
# SECTION 0: INSTALL & LOAD PACKAGES
# ============================================================================
# Run this block only once to install required packages:
# install.packages(c("dplyr", "plyr", "moments", "PerformanceAnalytics",
#                    "stabs", "ggplot2", "ggeffects", "readxl",
#                    "RColorBrewer", "corrplot", "gridExtra", "lme4"))

library(dplyr)
library(plyr)          # join()
library(moments)       # skewness()
library(PerformanceAnalytics)  # chart.Correlation()
library(stabs)         # stabsel() — stability selection
library(ggplot2)
library(ggeffects)     # ggpredict() for marginal effects plots
library(readxl)
library(RColorBrewer)
library(corrplot)
library(gridExtra)

# Set working directory — adjust to your folder path:
# setwd("C:/Users/YourName/Documents/LeafDecomp_2025")

# ============================================================================
# SECTION 1: LOAD ALL RAW DATA FILES
# ============================================================================

# 1a. Processed sample weights (AFDM measurements)
processed <- read_excel("2025_field_processed_samples.xlsx")
cat("Processed samples loaded:", nrow(processed), "rows\n")

# 1b. Initial bag masses before deployment
initial_mass <- read_excel("2025_field_initial_mass.xlsx")
cat("Initial masses loaded:", nrow(initial_mass), "rows\n")

# 1c. Transport controls (handling loss)
transport_ctrl <- read_excel("2025_field_transport_controls.xlsx")
cat("Transport controls loaded:", nrow(transport_ctrl), "rows\n")

# 1d. Combined environmental variables (water chemistry, substrate, land use, invertebrates)
env_vars <- read_excel("SFB_Combined_Data_Test.xlsx", sheet = "All_Data")
cat("Environmental variables loaded:", nrow(env_vars), "rows,", ncol(env_vars), "columns\n")

# 1e. Temperature logger data (Boye)
temp_logger <- read.csv("monthly_Temp_logger_BOYE_march2025.csv")
cat("Temperature logger loaded:", nrow(temp_logger), "rows\n")

# 1f. Field date summary (deployment and retrieval dates)
field_dates <- read_excel("Boye_Kinzig_Field_Date_Summary.xlsx", sheet = "Bag_Details_All")
cat("Field dates loaded:", nrow(field_dates), "rows\n")

# ============================================================================
# SECTION 2: COMPUTE AFDM END (Ash-Free Dry Mass at retrieval)
# ============================================================================
# AFDM = mass dried at 60°C  −  mass combusted at 500°C
# This gives the organic fraction only (ash = inorganic mineral fraction)

# Disk correction: account for material removed as leaf disks for other analyses
# Standard disk weight = 12 mg = 0.012 g per disk
DISK_WEIGHT_G <- 0.012

processed <- processed %>%
  mutate(
    afdm_end = rest_60_g - rest_500_g,
    no_disks_all = replace_na(no_disks_all, 0),
    afdm_end_corrected = afdm_end + no_disks_all * DISK_WEIGHT_G
  )

cat("AFDM end computed. Summary:\n")
print(summary(processed$afdm_end_corrected))

# ============================================================================
# SECTION 3: COMPUTE AFDM START (Initial AFDM from transport controls)
# ============================================================================
# Transport controls quantify handling/transport losses.
# Use them to convert air-dry initial mass to theoretical AFDM start.

transport_ctrl <- transport_ctrl %>%
  mutate(afdm_transport = mass_60_g - mass_500_g)

# Compute shrinkage ratios per catchment
get_ratios <- function(tc) {
  tc_clean <- tc %>% filter(!is.na(mass_60_g), !is.na(mass_500_g), !is.na(initial_mass_g))
  perc_60  <- mean(tc_clean$mass_60_g  / tc_clean$initial_mass_g)
  perc_500 <- mean(tc_clean$mass_500_g / tc_clean$initial_mass_g)
  return(list(perc_60 = perc_60, perc_500 = perc_500))
}

ratios_boye   <- get_ratios(filter(transport_ctrl, catchment == "Boye"))
ratios_kinzig <- get_ratios(filter(transport_ctrl, catchment == "Kinzig"))

cat(sprintf("Boye ratios: perc_60 = %.4f, perc_500 = %.4f\n",
            ratios_boye$perc_60, ratios_boye$perc_500))
cat(sprintf("Kinzig ratios: perc_60 = %.4f, perc_500 = %.4f\n",
            ratios_kinzig$perc_60, ratios_kinzig$perc_500))

# Merge initial mass with processed data
processed <- processed %>%
  left_join(initial_mass %>% select(bag_ID, initial_mass_g), by = "bag_ID") %>%
  mutate(
    afdm_start = case_when(
      catchment == "Boye"   ~ initial_mass_g * (ratios_boye$perc_60   - ratios_boye$perc_500),
      catchment == "Kinzig" ~ initial_mass_g * (ratios_kinzig$perc_60 - ratios_kinzig$perc_500),
      TRUE ~ NA_real_
    )
  )

cat("AFDM start computed. Summary:\n")
print(summary(processed$afdm_start))

# ============================================================================
# SECTION 4: COMPUTE INCUBATION DAYS AND DEGREE-DAYS
# ============================================================================

processed <- processed %>%
  mutate(
    date_deployment = as.Date(date_deployment),
    date_sampled    = as.Date(date_sampled),
    incubation_days = as.numeric(date_sampled - date_deployment)
  )

# Merge temperature data
# From temperature logger (Boye) and environmental data (Kinzig + remaining Boye sites)
site_temp <- env_vars %>%
  select(site, avg_temp_degC) %>%
  filter(!is.na(avg_temp_degC)) %>%
  distinct(site, .keep_all = TRUE)

# Also use the logger data
temp_logger_df <- temp_logger %>%
  rename(stream = Stream, avg_temp_logger = Average) %>%
  select(stream, avg_temp_logger)

processed <- processed %>%
  left_join(site_temp, by = c("stream" = "site")) %>%
  left_join(temp_logger_df, by = "stream") %>%
  mutate(
    temp = coalesce(avg_temp_degC, avg_temp_logger,
                    ifelse(catchment == "Boye", 8.0, 10.0)),
    degree_days = temp * incubation_days
  )

cat("Incubation days summary:\n")
print(summary(processed$incubation_days))
cat("Degree-days summary:\n")
print(summary(processed$degree_days))

# ============================================================================
# SECTION 5: COMPUTE DSM (Decomposed Substrate Mass)
# ============================================================================
# DSM formula (Schreiner et al. 2023):
# DSM = (AFDM_start - AFDM_end) / AFDM_start / degree_days × 100
#
# This gives % of initial organic mass decomposed per degree-day.
# Negative values indicate no decomposition and are removed.

processed <- processed %>%
  mutate(
    mass_loss_g = afdm_start - afdm_end_corrected,
    DSM = (mass_loss_g / afdm_start) / degree_days * 100,
    DSM = ifelse(DSM <= 0, NA_real_, DSM)  # Remove non-positive DSM
  )

cat("\nDSM computed. Summary:\n")
print(processed %>%
  group_by(catchment, type) %>%
  summarise(
    n      = sum(!is.na(DSM)),
    mean   = round(mean(DSM, na.rm = TRUE), 5),
    sd     = round(sd(DSM,   na.rm = TRUE), 5),
    min    = round(min(DSM,   na.rm = TRUE), 5),
    max    = round(max(DSM,   na.rm = TRUE), 5),
    .groups = "drop"
  ))

# Save processed data
write.csv(processed, "LeafDecomp_Output/data/processed_DSM_data.csv", row.names = FALSE)

# ============================================================================
# SECTION 6: SITE-LEVEL DSM MEANS AND ENVIRONMENTAL VARIABLES
# ============================================================================

# Compute mean DSM per site and bag type
dsm_site <- processed %>%
  filter(!is.na(DSM)) %>%
  group_by(stream, catchment, type) %>%
  summarise(
    DSM_mean = mean(DSM, na.rm = TRUE),
    DSM_sd   = sd(DSM,   na.rm = TRUE),
    DSM_n    = n(),
    .groups  = "drop"
  )

# Pivot to wide format
dsm_wide <- dsm_site %>%
  tidyr::pivot_wider(
    id_cols     = c(stream, catchment),
    names_from  = type,
    values_from = c(DSM_mean, DSM_sd, DSM_n)
  )

# Merge environmental variables
env_select <- env_vars %>%
  select(site, pH_dimensionless, conductivity_uS_cm, dissolved_oxygen_mg_l,
         oxygen_saturation_pct, water_temp_degC, water_level_cm,
         ortho_phosphate_mg_l, total_phosphate_mg_l, ammonia_mg_l,
         nitrite_mg_l, nitrate_mg_l, total_nitrogen_mg_l,
         chloride_mg_l, sulfate_mg_l,
         Cd_mg_l, Cu_mg_l, Fe_mg_l, Ni_mg_l, Pb_mg_l, Zn_mg_l,
         Ca_mg_l, Mg_mg_l, Na_mg_l, K_mg_l, Co3_mg_l, HCO3_mg_l,
         avg_temp_degC, sd_temp,
         makro, meso, mikro, akal, psam, tech_1, submers, emers,
         `living parts`, xylal, CPOM, FPOM,
         total_area, pct_agricultural, pct_industrial, pct_urban, pct_rural,
         Gammarus_pulex, Gammarus_fossarum, `Gammarus sp.`)

# Combine Gammarus species
env_select <- env_select %>%
  mutate(Gammarus_total = rowSums(select(., starts_with("Gammarus")), na.rm = TRUE))

master_data <- dsm_wide %>%
  left_join(env_select, by = c("stream" = "site"))

cat("\nMaster dataset: ", nrow(master_data), "sites,", ncol(master_data), "variables\n")

# Split by stream
boye_data   <- filter(master_data, catchment == "Boye")
kinzig_data <- filter(master_data, catchment == "Kinzig")
cat("Boye: ", nrow(boye_data), "sites\n")
cat("Kinzig: ", nrow(kinzig_data), "sites\n")

# ============================================================================
# SECTION 7: EXPLORATORY DATA ANALYSIS
# ============================================================================

# ---- 7a. Check skewness of key variables (Boye)
cat("\n=== Variable Skewness Check (Boye) ===\n")
boye_numeric <- boye_data %>%
  select(DSM_mean_fine, DSM_mean_coarse, conductivity_uS_cm, nitrate_mg_l,
         total_nitrogen_mg_l, Zn_mg_l, Cu_mg_l, Ni_mg_l, pct_urban, Gammarus_total) %>%
  na.omit()

skew_results <- sapply(boye_numeric, skewness)
print(round(skew_results, 3))
cat("Variables with |skewness| > 1 should be log-transformed\n")

# ---- 7b. Log-transform skewed variables
# Log1p = log(x + 1) — safe for zeros
log_vars <- c("Zn_mg_l","Cu_mg_l","Ni_mg_l","Pb_mg_l","Cd_mg_l",
              "Gammarus_total","total_invertebrate_abundance",
              "ortho_phosphate_mg_l","total_phosphate_mg_l","ammonia_mg_l")

for (col in log_vars) {
  if (col %in% names(boye_data)) {
    boye_data[[paste0("log_", col)]] <- log1p(boye_data[[col]])
  }
  if (col %in% names(kinzig_data)) {
    kinzig_data[[paste0("log_", col)]] <- log1p(kinzig_data[[col]])
  }
}

# ---- 7c. Correlation matrix (Boye)
cat("\n=== Pearson Correlation Matrix — Boye (key variables) ===\n")
boye_corr_vars <- boye_data %>%
  select(DSM_mean_fine, DSM_mean_coarse, conductivity_uS_cm, pH_dimensionless,
         nitrate_mg_l, total_nitrogen_mg_l, ortho_phosphate_mg_l,
         log_Zn_mg_l, log_Cu_mg_l, Ni_mg_l, pct_urban, pct_agricultural,
         log_Gammarus_total, avg_temp_degC) %>%
  na.omit()

print(round(cor(boye_corr_vars, use = "complete.obs"), 3))

# ---- 7d. VIF check to identify collinear variables (Boye fine model)
library(car)
# Full model to check VIF (must have ≥3 observations per variable)
vif_model_boye <- lm(DSM_mean_fine ~ conductivity_uS_cm + pH_dimensionless +
                       nitrate_mg_l + log_Zn_mg_l + pct_urban + avg_temp_degC,
                     data = boye_data)
cat("\nVIF for Boye fine DSM model:\n")
print(vif(vif_model_boye))
cat("Remove variables with VIF > 10 from models\n")

# ============================================================================
# SECTION 8: HYPOTHESIS TESTING — BOYE STREAM
# ============================================================================
# Four hypotheses for Boye (urban/industrial catchment):
#
# H1: Higher conductivity is negatively associated with fine mesh DSM
# H2: Higher urban land use percentage reduces leaf decomposition rates
# H3: Higher Gammarus shredder abundance drives higher coarse mesh DSM
# H4: Higher Zn and Cu concentrations reduce fine mesh (microbial) DSM

cat("\n========== BOYE STREAM — HYPOTHESIS TESTING ==========\n")

# ----- H1: Conductivity vs DSM fine -----
cat("\n--- H1: Conductivity vs Fine Mesh DSM ---\n")
h1_boye <- cor.test(boye_data$conductivity_uS_cm, boye_data$DSM_mean_fine,
                    use = "complete.obs", method = "pearson")
cat("Pearson r:", round(h1_boye$estimate, 3), "\n")
cat("p-value:", round(h1_boye$p.value, 4), "\n")
cat("Interpretation:", ifelse(h1_boye$p.value < 0.05,
    "SIGNIFICANT — conductivity is associated with fine DSM",
    "NOT SIGNIFICANT"), "\n")

lm_h1_boye <- lm(DSM_mean_fine ~ conductivity_uS_cm, data = boye_data)
print(summary(lm_h1_boye))

# ----- H2: Urban land use vs DSM -----
cat("\n--- H2: Urban Land Use vs Decomposition Rates ---\n")
h2_fine <- cor.test(boye_data$pct_urban, boye_data$DSM_mean_fine,
                    use = "complete.obs", method = "pearson")
h2_coarse <- cor.test(boye_data$pct_urban, boye_data$DSM_mean_coarse,
                      use = "complete.obs", method = "pearson")
cat("Fine mesh — r:", round(h2_fine$estimate,3), "p:", round(h2_fine$p.value,4), "\n")
cat("Coarse mesh — r:", round(h2_coarse$estimate,3), "p:", round(h2_coarse$p.value,4), "\n")

# Linear models
lm_h2_fine   <- lm(DSM_mean_fine   ~ pct_urban, data = boye_data)
lm_h2_coarse <- lm(DSM_mean_coarse ~ pct_urban, data = boye_data)
cat("Fine mesh R²:", round(summary(lm_h2_fine)$r.squared, 3), "\n")
cat("Coarse mesh R²:", round(summary(lm_h2_coarse)$r.squared, 3), "\n")

# ----- H3: Gammarus abundance vs coarse DSM -----
cat("\n--- H3: Gammarus Shredder Abundance vs Coarse DSM ---\n")
boye_h3 <- boye_data %>% filter(!is.na(Gammarus_total), !is.na(DSM_mean_coarse))
h3_cor <- cor.test(log1p(boye_h3$Gammarus_total), boye_h3$DSM_mean_coarse,
                   method = "pearson")
cat("Log-Gammarus vs Coarse DSM — r:", round(h3_cor$estimate,3),
    "p:", round(h3_cor$p.value,4), "\n")

lm_h3 <- lm(DSM_mean_coarse ~ log1p(Gammarus_total), data = boye_h3)
print(summary(lm_h3))

# Compare to fine mesh (H3 prediction: coarse > fine due to shredders)
t_test_boye <- t.test(
  filter(processed, catchment == "Boye", type == "coarse")$DSM,
  filter(processed, catchment == "Boye", type == "fine")$DSM,
  na.action = na.omit
)
cat("\nCoarse vs Fine DSM (Boye) t-test:\n")
cat("t =", round(t_test_boye$statistic, 3), "\n")
cat("p =", round(t_test_boye$p.value, 6), "\n")
cat("Mean Coarse:", round(t_test_boye$estimate[1], 5),
    " | Mean Fine:", round(t_test_boye$estimate[2], 5), "\n")

# ----- H4: Metal contamination vs fine DSM -----
cat("\n--- H4: Metal Contamination (Zn, Cu) vs Fine Mesh DSM ---\n")
h4_zn <- cor.test(boye_data$Zn_mg_l, boye_data$DSM_mean_fine,
                  use = "complete.obs", method = "pearson")
h4_cu <- cor.test(boye_data$Cu_mg_l, boye_data$DSM_mean_fine,
                  use = "complete.obs", method = "pearson")
cat("Zinc vs Fine DSM — r:", round(h4_zn$estimate,3), "p:", round(h4_zn$p.value,4), "\n")
cat("Copper vs Fine DSM — r:", round(h4_cu$estimate,3), "p:", round(h4_cu$p.value,4), "\n")

# Multiple regression model
lm_h4 <- lm(DSM_mean_fine ~ log1p(Zn_mg_l) + log1p(Cu_mg_l), data = boye_data)
print(summary(lm_h4))

# ============================================================================
# SECTION 9: HYPOTHESIS TESTING — KINZIG STREAM
# ============================================================================
# Four hypotheses for Kinzig (agricultural/rural catchment):
#
# H1: Higher nitrate concentrations (from agricultural runoff) are positively
#     correlated with fine mesh DSM (nutrient-stimulated microbial activity)
# H2: Higher agricultural land use percentage is positively associated with
#     DSM due to nutrient-enhanced microbial decomposition
# H3: Coarse mesh DSM is significantly higher than fine mesh DSM,
#     indicating invertebrate shredder contribution to decomposition
# H4: Higher conductivity (nutrient-enriched sites) is positively
#     associated with higher leaf decomposition rates

cat("\n========== KINZIG STREAM — HYPOTHESIS TESTING ==========\n")

# ----- H1: Nitrate vs fine DSM -----
cat("\n--- H1: Nitrate vs Fine Mesh DSM ---\n")
h1_kinzig <- cor.test(kinzig_data$nitrate_mg_l, kinzig_data$DSM_mean_fine,
                      use = "complete.obs", method = "pearson")
cat("Pearson r:", round(h1_kinzig$estimate, 3), "\n")
cat("p-value:", round(h1_kinzig$p.value, 4), "\n")
cat("Interpretation:", ifelse(h1_kinzig$p.value < 0.05,
    "SIGNIFICANT — nitrate associated with fine mesh DSM",
    "NOT SIGNIFICANT"), "\n")

lm_h1_kin <- lm(DSM_mean_fine ~ nitrate_mg_l, data = kinzig_data)
print(summary(lm_h1_kin))

# ----- H2: Agricultural land use vs DSM -----
cat("\n--- H2: Agricultural Land Use vs Decomposition Rates ---\n")
h2k_fine   <- cor.test(kinzig_data$pct_agricultural, kinzig_data$DSM_mean_fine,
                        use = "complete.obs", method = "pearson")
h2k_coarse <- cor.test(kinzig_data$pct_agricultural, kinzig_data$DSM_mean_coarse,
                        use = "complete.obs", method = "pearson")
cat("Fine mesh — r:", round(h2k_fine$estimate,3), "p:", round(h2k_fine$p.value,4), "\n")
cat("Coarse mesh — r:", round(h2k_coarse$estimate,3), "p:", round(h2k_coarse$p.value,4), "\n")

lm_h2k <- lm(DSM_mean_fine ~ pct_agricultural, data = kinzig_data)
cat("Agricultural land use fine DSM R²:", round(summary(lm_h2k)$r.squared, 3), "\n")

# ----- H3: Fine vs Coarse DSM comparison -----
cat("\n--- H3: Fine vs Coarse Mesh DSM (Shredder Contribution) ---\n")
t_test_kin <- t.test(
  filter(processed, catchment == "Kinzig", type == "coarse")$DSM,
  filter(processed, catchment == "Kinzig", type == "fine")$DSM,
  na.action = na.omit
)
cat("t =", round(t_test_kin$statistic, 3), "\n")
cat("p =", round(t_test_kin$p.value, 6), "\n")
cat("Mean Coarse:", round(t_test_kin$estimate[1], 5),
    " | Mean Fine:", round(t_test_kin$estimate[2], 5), "\n")
cat("Shredder effect (Coarse - Fine):",
    round(t_test_kin$estimate[1] - t_test_kin$estimate[2], 5), "\n")

# ----- H4: Conductivity vs DSM -----
cat("\n--- H4: Conductivity vs Decomposition Rates ---\n")
h4k_fine   <- cor.test(kinzig_data$conductivity_uS_cm, kinzig_data$DSM_mean_fine,
                        use = "complete.obs", method = "pearson")
h4k_coarse <- cor.test(kinzig_data$conductivity_uS_cm, kinzig_data$DSM_mean_coarse,
                        use = "complete.obs", method = "pearson")
cat("Fine mesh — r:", round(h4k_fine$estimate,3), "p:", round(h4k_fine$p.value,4), "\n")
cat("Coarse mesh — r:", round(h4k_coarse$estimate,3), "p:", round(h4k_coarse$p.value,4), "\n")

# ============================================================================
# SECTION 10: MULTIPLE REGRESSION MODELS
# ============================================================================

cat("\n========== MULTIPLE REGRESSION MODELS ==========\n")

# ---- 10a. Boye Fine Mesh — key predictors ----
cat("\n--- Boye Fine Mesh DSM — Multiple Regression ---\n")
boye_model_fine <- lm(DSM_mean_fine ~ conductivity_uS_cm + log1p(Zn_mg_l) +
                        pct_urban + avg_temp_degC,
                      data = boye_data)
print(summary(boye_model_fine))

# ---- 10b. Boye Coarse Mesh ----
cat("\n--- Boye Coarse Mesh DSM — Multiple Regression ---\n")
boye_model_coarse <- lm(DSM_mean_coarse ~ log1p(Gammarus_total) +
                           conductivity_uS_cm + pct_urban,
                         data = boye_data)
print(summary(boye_model_coarse))

# ---- 10c. Kinzig Fine Mesh ----
cat("\n--- Kinzig Fine Mesh DSM — Multiple Regression ---\n")
kinzig_model_fine <- lm(DSM_mean_fine ~ nitrate_mg_l + conductivity_uS_cm +
                           pct_agricultural + avg_temp_degC,
                         data = kinzig_data)
print(summary(kinzig_model_fine))

# ---- 10d. Kinzig Coarse Mesh ----
cat("\n--- Kinzig Coarse Mesh DSM — Multiple Regression ---\n")
kinzig_model_coarse <- lm(DSM_mean_coarse ~ nitrate_mg_l + conductivity_uS_cm +
                              pct_agricultural,
                            data = kinzig_data)
print(summary(kinzig_model_coarse))

# ============================================================================
# SECTION 11: STABILITY SELECTION (stabs package)
# ============================================================================
# Stability selection identifies robust predictors across resampled sub-models.
# Variables selected in > 50% of runs are considered stable predictors.
# This is the approach used in Schreiner et al. (2023).

cat("\n========== STABILITY SELECTION ==========\n")

# --- Prepare predictor matrices ---
prep_pred_matrix <- function(df, pred_cols) {
  df_sub <- df %>% select(all_of(pred_cols)) %>% na.omit()
  scale(df_sub)  # standardize all predictors to z-scores
}

# Boye Fine predictors (remove highly collinear variables, VIF > 10)
boye_pred_cols_fine <- c("conductivity_uS_cm", "pH_dimensionless",
                          "nitrate_mg_l", "log_Zn_mg_l", "log_Cu_mg_l",
                          "pct_urban", "avg_temp_degC")

# Only run if stabs is available and data has sufficient rows
if (requireNamespace("stabs", quietly = TRUE) && nrow(boye_data) > 10) {
  boye_pred_fine <- prep_pred_matrix(
    boye_data %>% mutate(log_Zn_mg_l = log1p(Zn_mg_l), log_Cu_mg_l = log1p(Cu_mg_l)),
    boye_pred_cols_fine
  )
  # Match response to predictor rows
  complete_idx_bf <- which(complete.cases(
    boye_data %>% select(DSM_mean_fine, all_of(boye_pred_cols_fine))
  ))
  boye_dsm_fine_complete <- boye_data$DSM_mean_fine[complete_idx_bf]

  set.seed(42)
  stabsel_boye_fine <- stabsel(
    x     = boye_pred_fine,
    y     = boye_dsm_fine_complete,
    fitfun = glmnet.lasso,
    cutoff = 0.6,
    PFER   = 1
  )
  cat("\nBoye Fine DSM — Stability Selected Variables:\n")
  print(stabsel_boye_fine$selected)
  print(stabsel_boye_fine$max)
}

# Kinzig Fine predictors
kinzig_pred_cols_fine <- c("conductivity_uS_cm", "pH_dimensionless",
                            "nitrate_mg_l", "total_nitrogen_mg_l",
                            "pct_agricultural", "avg_temp_degC")

if (requireNamespace("stabs", quietly = TRUE) && nrow(kinzig_data) > 10) {
  kinzig_pred_fine <- prep_pred_matrix(kinzig_data, kinzig_pred_cols_fine)
  complete_idx_kf <- which(complete.cases(
    kinzig_data %>% select(DSM_mean_fine, all_of(kinzig_pred_cols_fine))
  ))
  kinzig_dsm_fine_complete <- kinzig_data$DSM_mean_fine[complete_idx_kf]

  set.seed(42)
  stabsel_kinzig_fine <- stabsel(
    x     = kinzig_pred_fine,
    y     = kinzig_dsm_fine_complete,
    fitfun = glmnet.lasso,
    cutoff = 0.6,
    PFER   = 1
  )
  cat("\nKinzig Fine DSM — Stability Selected Variables:\n")
  print(stabsel_kinzig_fine$selected)
  print(stabsel_kinzig_fine$max)
}

# ============================================================================
# SECTION 12: VISUALIZATION — PUBLICATION-QUALITY PLOTS
# ============================================================================

cat("\n========== GENERATING PLOTS ==========\n")
plot_dir <- "LeafDecomp_Output/plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- PLOT 1: DSM comparison — fine vs coarse, both streams ----
# Prepare boxplot data
plot_data <- processed %>%
  filter(!is.na(DSM)) %>%
  mutate(
    stream_type = paste(catchment, type, sep = " — "),
    type = factor(type, levels = c("fine","coarse"), labels = c("Fine mesh","Coarse mesh")),
    catchment = factor(catchment, levels = c("Boye","Kinzig"))
  )

p1 <- ggplot(plot_data, aes(x = type, y = DSM, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5, outlier.alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  facet_wrap(~ catchment, scales = "free_y") +
  scale_fill_manual(values = c("Fine mesh" = "#3498db", "Coarse mesh" = "#e74c3c")) +
  labs(
    title    = "Leaf Decomposition Rates (DSM) by Stream and Bag Type",
    subtitle = "DSM = % organic mass lost per degree-day",
    x        = "Bag type",
    y        = "DSM (% mass loss · degree-day⁻¹)",
    fill     = "Bag type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#2E75B6", color = "#2E75B6"),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13)
  )
ggsave(file.path(plot_dir, "R_Fig1_DSM_FineVsCoarse_Boxplot.pdf"), p1,
       width = 10, height = 6, dpi = 300)
cat("R_Fig1 saved\n")

# ---- PLOT 2: DSM per site — Boye ----
boye_site_dsm <- processed %>%
  filter(catchment == "Boye", !is.na(DSM)) %>%
  group_by(stream, type) %>%
  summarise(mean_DSM = mean(DSM, na.rm=TRUE),
            sd_DSM   = sd(DSM, na.rm=TRUE), .groups="drop") %>%
  mutate(type = factor(type, levels=c("fine","coarse"), labels=c("Fine","Coarse")))

p2 <- ggplot(boye_site_dsm, aes(x = mean_DSM, y = reorder(stream, mean_DSM), fill = type)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_errorbarh(aes(xmin = pmax(mean_DSM - sd_DSM, 0), xmax = mean_DSM + sd_DSM),
                 position = position_dodge(0.9), height = 0.3) +
  scale_fill_manual(values = c("Fine" = "#3498db", "Coarse" = "#e74c3c")) +
  labs(title = "Boye Stream — DSM per Site",
       x = "Mean DSM (% loss · degree-day⁻¹)", y = "Site", fill = "Bag type") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face="bold"))
ggsave(file.path(plot_dir, "R_Fig2_Boye_DSM_perSite.pdf"), p2,
       width = 10, height = 8, dpi = 300)
cat("R_Fig2 saved\n")

# ---- PLOT 3: DSM per site — Kinzig ----
kinzig_site_dsm <- processed %>%
  filter(catchment == "Kinzig", !is.na(DSM)) %>%
  group_by(stream, type) %>%
  summarise(mean_DSM = mean(DSM, na.rm=TRUE),
            sd_DSM   = sd(DSM, na.rm=TRUE), .groups="drop") %>%
  mutate(type = factor(type, levels=c("fine","coarse"), labels=c("Fine","Coarse")))

p3 <- ggplot(kinzig_site_dsm, aes(x = mean_DSM, y = reorder(stream, mean_DSM), fill = type)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_errorbarh(aes(xmin = pmax(mean_DSM - sd_DSM, 0), xmax = mean_DSM + sd_DSM),
                 position = position_dodge(0.9), height = 0.3) +
  scale_fill_manual(values = c("Fine" = "#27ae60", "Coarse" = "#f39c12")) +
  labs(title = "Kinzig Stream — DSM per Site",
       x = "Mean DSM (% loss · degree-day⁻¹)", y = "Site", fill = "Bag type") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face="bold"))
ggsave(file.path(plot_dir, "R_Fig3_Kinzig_DSM_perSite.pdf"), p3,
       width = 10, height = 8, dpi = 300)
cat("R_Fig3 saved\n")

# ---- PLOT 4: Scatter regression — Boye H1 (Conductivity vs DSM) ----
p4 <- ggplot(boye_data %>% filter(!is.na(conductivity_uS_cm), !is.na(DSM_mean_fine)),
             aes(x = conductivity_uS_cm, y = DSM_mean_fine)) +
  geom_point(color = "#2E75B6", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80", alpha = 0.3) +
  labs(title = "Boye Stream — H1: Conductivity vs Fine Mesh DSM",
       subtitle = paste("Pearson r =", round(h1_boye$estimate,3),
                        "| p =", round(h1_boye$p.value, 4)),
       x = "Conductivity (µS/cm)", y = "Mean DSM — Fine mesh") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"))
ggsave(file.path(plot_dir, "R_Fig4_Boye_H1_Conductivity.pdf"), p4,
       width = 7, height = 5, dpi = 300)
cat("R_Fig4 saved\n")

# ---- PLOT 5: Boye H3 — Gammarus vs Coarse DSM ----
p5 <- ggplot(boye_data %>% filter(!is.na(Gammarus_total), !is.na(DSM_mean_coarse)),
             aes(x = Gammarus_total + 1, y = DSM_mean_coarse)) +
  geom_point(color = "#e74c3c", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE,
              color = "black", fill = "grey80", alpha = 0.3) +
  scale_x_log10(labels = scales::comma) +
  labs(title = "Boye Stream — H3: Gammarus Abundance vs Coarse DSM",
       subtitle = paste("Pearson r (log) =", round(h3_cor$estimate,3),
                        "| p =", round(h3_cor$p.value, 4)),
       x = "Gammarus abundance (ind/m², log scale)", y = "Mean DSM — Coarse mesh") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"))
ggsave(file.path(plot_dir, "R_Fig5_Boye_H3_Gammarus.pdf"), p5,
       width = 7, height = 5, dpi = 300)
cat("R_Fig5 saved\n")

# ---- PLOT 6: Kinzig H1 — Nitrate vs Fine DSM ----
p6 <- ggplot(kinzig_data %>% filter(!is.na(nitrate_mg_l), !is.na(DSM_mean_fine)),
             aes(x = nitrate_mg_l, y = DSM_mean_fine)) +
  geom_point(color = "#27ae60", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80", alpha = 0.3) +
  labs(title = "Kinzig Stream — H1: Nitrate vs Fine Mesh DSM",
       subtitle = paste("Pearson r =", round(h1_kinzig$estimate,3),
                        "| p =", round(h1_kinzig$p.value, 4)),
       x = "Nitrate (mg/L)", y = "Mean DSM — Fine mesh") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold"))
ggsave(file.path(plot_dir, "R_Fig6_Kinzig_H1_Nitrate.pdf"), p6,
       width = 7, height = 5, dpi = 300)
cat("R_Fig6 saved\n")

# ---- PLOT 7: Kinzig H3 — Fine vs Coarse by site ----
p7 <- ggplot(kinzig_site_dsm %>% filter(!is.na(mean_DSM)),
             aes(x = stream, y = mean_DSM, fill = type)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(mean_DSM - sd_DSM, 0), ymax = mean_DSM + sd_DSM),
                position = position_dodge(0.9), width = 0.25) +
  scale_fill_manual(values = c("Fine" = "#27ae60", "Coarse" = "#f39c12")) +
  labs(title = "Kinzig Stream — H3: Fine vs Coarse DSM (Shredder Effect)",
       subtitle = paste("Coarse > Fine: t-test p =", round(t_test_kin$p.value, 6)),
       x = "Site", y = "Mean DSM", fill = "Bag type") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face="bold"), legend.position = "bottom")
ggsave(file.path(plot_dir, "R_Fig7_Kinzig_H3_FineVsCoarse_Sites.pdf"), p7,
       width = 11, height = 6, dpi = 300)
cat("R_Fig7 saved\n")

# ---- PLOT 8: Marginal effects from Boye fine model ----
if (exists("boye_model_fine")) {
  pred_plots_boye <- ggpredict(boye_model_fine, terms = c("conductivity_uS_cm"))
  p8 <- plot(pred_plots_boye) +
    labs(title = "Marginal Effect: Conductivity on Boye Fine Mesh DSM",
         x = "Conductivity (µS/cm)", y = "Predicted DSM (fine mesh)") +
    theme_bw(base_size = 12)
  ggsave(file.path(plot_dir, "R_Fig8_Boye_MarginalEffect_Conductivity.pdf"), p8,
         width = 7, height = 5, dpi = 300)
  cat("R_Fig8 saved\n")
}

# ---- PLOT 9: Correlation heatmap — Boye ----
pdf(file.path(plot_dir, "R_Fig9_Boye_Correlation_Heatmap.pdf"), width=12, height=10)
boye_corr_mat <- cor(boye_corr_vars, use = "complete.obs")
corrplot(boye_corr_mat, method = "color", type = "lower",
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         col = colorRampPalette(c("#c0392b","white","#2E75B6"))(200),
         title = "Boye Stream — Pearson Correlation Matrix",
         mar = c(0,0,2,0))
dev.off()
cat("R_Fig9 saved\n")

# ---- PLOT 10: Summary overview ----
summary_df <- processed %>%
  filter(!is.na(DSM)) %>%
  group_by(catchment, type) %>%
  summarise(mean = mean(DSM), sd = sd(DSM), .groups="drop") %>%
  mutate(label = paste(catchment, "\n", type))

p10 <- ggplot(summary_df, aes(x = label, y = mean, fill = interaction(catchment, type))) +
  geom_col(alpha = 0.85, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(mean-sd, 0), ymax = mean+sd), width = 0.2) +
  geom_text(aes(label = round(mean, 3)), vjust = -0.5, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("Boye.fine"="#3498db","Boye.coarse"="#e74c3c",
                                "Kinzig.fine"="#27ae60","Kinzig.coarse"="#f39c12")) +
  labs(title = "Summary: Mean DSM by Stream and Bag Type",
       x = "", y = "Mean DSM ± SD (% loss · degree-day⁻¹)", fill = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face="bold"))
ggsave(file.path(plot_dir, "R_Fig10_DSM_Summary_Overview.pdf"), p10,
       width = 8, height = 5, dpi = 300)
cat("R_Fig10 saved\n")

# ============================================================================
# SECTION 13: EXPORT RESULTS TABLES
# ============================================================================

cat("\n========== EXPORTING RESULTS ==========\n")

# Summary table
results_summary <- processed %>%
  filter(!is.na(DSM)) %>%
  group_by(catchment, type) %>%
  summarise(
    n_bags      = n(),
    n_sites     = n_distinct(stream),
    mean_DSM    = round(mean(DSM), 6),
    sd_DSM      = round(sd(DSM), 6),
    median_DSM  = round(median(DSM), 6),
    min_DSM     = round(min(DSM), 6),
    max_DSM     = round(max(DSM), 6),
    .groups = "drop"
  )

write.csv(results_summary, "LeafDecomp_Output/results/DSM_summary_statistics.csv",
          row.names = FALSE)
print(results_summary)

# Site-level results
write.csv(master_data, "LeafDecomp_Output/results/master_site_level_data.csv",
          row.names = FALSE)
cat("All results exported to LeafDecomp_Output/results/\n")

# ============================================================================
# SECTION 14: HYPOTHESIS RESULTS SUMMARY TABLE
# ============================================================================

hypothesis_results <- data.frame(
  Stream = c(rep("Boye", 8), rep("Kinzig", 8)),
  Hypothesis = c(
    "H1","H1","H2","H2","H3","H3","H4","H4",
    "H1","H1","H2","H2","H3","H3","H4","H4"
  ),
  Predictor = c(
    "Conductivity","Conductivity","% Urban","% Urban",
    "Gammarus abundance","Gammarus abundance","Zinc (Zn)","Copper (Cu)",
    "Nitrate","Total Nitrogen","% Agricultural","% Agricultural",
    "Coarse vs Fine","Coarse vs Fine","Conductivity","Conductivity"
  ),
  Response = c(
    "Fine DSM","Coarse DSM","Fine DSM","Coarse DSM",
    "Coarse DSM","Fine DSM (control)","Fine DSM","Fine DSM",
    "Fine DSM","Fine DSM","Fine DSM","Coarse DSM",
    "Coarse DSM","Fine DSM","Fine DSM","Coarse DSM"
  ),
  Pearson_r = c(
    round(h1_boye$estimate,3), NA, round(h2_fine$estimate,3), round(h2_coarse$estimate,3),
    round(h3_cor$estimate,3), NA, round(h4_zn$estimate,3), round(h4_cu$estimate,3),
    round(h1_kinzig$estimate,3), NA, round(h2k_fine$estimate,3), round(h2k_coarse$estimate,3),
    NA, NA, round(h4k_fine$estimate,3), round(h4k_coarse$estimate,3)
  ),
  p_value = c(
    round(h1_boye$p.value,4), NA, round(h2_fine$p.value,4), round(h2_coarse$p.value,4),
    round(h3_cor$p.value,4), NA, round(h4_zn$p.value,4), round(h4_cu$p.value,4),
    round(h1_kinzig$p.value,4), NA, round(h2k_fine$p.value,4), round(h2k_coarse$p.value,4),
    round(t_test_kin$p.value,6), NA, round(h4k_fine$p.value,4), round(h4k_coarse$p.value,4)
  ),
  Significant = ifelse(c(
    h1_boye$p.value, NA, h2_fine$p.value, h2_coarse$p.value,
    h3_cor$p.value, NA, h4_zn$p.value, h4_cu$p.value,
    h1_kinzig$p.value, NA, h2k_fine$p.value, h2k_coarse$p.value,
    t_test_kin$p.value, NA, h4k_fine$p.value, h4k_coarse$p.value
  ) < 0.05, "Yes *", "No", na="—"),
  Supported = c(
    ifelse(h1_boye$p.value < 0.05, "Partially","No"), "—",
    ifelse(h2_fine$p.value < 0.05, "Yes","No"), ifelse(h2_coarse$p.value < 0.05, "Yes","No"),
    ifelse(h3_cor$p.value < 0.05, "Yes","No"), "—",
    ifelse(h4_zn$p.value < 0.05, "Yes","No"), ifelse(h4_cu$p.value < 0.05, "Yes","No"),
    ifelse(h1_kinzig$p.value < 0.05, "Yes","No"), "—",
    ifelse(h2k_fine$p.value < 0.05, "Yes","No"), ifelse(h2k_coarse$p.value < 0.05, "Yes","No"),
    ifelse(t_test_kin$p.value < 0.05, "Yes","No"), "—",
    ifelse(h4k_fine$p.value < 0.05, "Yes","No"), ifelse(h4k_coarse$p.value < 0.05, "Yes","No")
  )
)

write.csv(hypothesis_results, "LeafDecomp_Output/results/hypothesis_test_results.csv",
          row.names = FALSE)
cat("\n=== HYPOTHESIS RESULTS SUMMARY ===\n")
print(hypothesis_results)

cat("\n\n========== ANALYSIS COMPLETE ==========\n")
cat("All outputs saved to: LeafDecomp_Output/\n")
cat("  Data:    LeafDecomp_Output/data/\n")
cat("  Plots:   LeafDecomp_Output/plots/\n")
cat("  Results: LeafDecomp_Output/results/\n")
