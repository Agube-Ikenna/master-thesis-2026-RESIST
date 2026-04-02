# ==================================================================
# Leaf Decomposition in the Kinzig Catchment – Full R Workflow
# Author:Agube Ikenna Stephen
# Last updated: 04.07.2025
# Description: Mixed model analysis, environmental predictors, figures
# ==================================================================




setwd("C:/Users/agube/OneDrive/Documentos/r_programming")
getwd()
library(readxl)

# Load all files
processed <- read_excel("2025_field_processed_samples.xlsx")
initial_mass <- read_excel("2025_field_initial_mass.xlsx")
transport <- read_excel("2025_field_transport_controls.xlsx")
site_info <- read_excel("2025_field_site_information_Kinzig.xlsx")

install.packages("readxl")
library(readxl)

library(readxl)

setwd("C:/Users/agube/OneDrive/Documentos/r_programming")

library(tidyverse)
library(readxl)

# Load all five field data files
processed     <- read_excel("2025_field_processed_samples.xlsx")
initial_mass  <- read_excel("2025_field_initial_mass.xlsx")
transport     <- read_excel("2025_field_transport_controls.xlsx")
site_info     <- read_excel("2025_field_site_information_Kinzig.xlsx")
inverts       <- read_excel("2025_field_invertebrates_coarse_leaf_bags.xlsx")  

# Clean invertebrate dataset if needed
inverts <- inverts %>%
  rename(bag_ID = bag_ID_column_name)  # replace with correct name if needed

# Join inverts to processed sample data (coarse mesh only)
coarse_with_inverts <- processed %>%
  filter(type == "coarse", catchment == "Kinzig") %>%
  left_join(inverts, by = c("stream", "replicate", "bag_ID"))

install.packages(c(
  "tidyverse",    # Core pipes, wrangling, plotting
  "readxl",       # Excel reading
  "lme4",         # Mixed models
  "lmerTest",     # p-values for lme4
  "emmeans",      # Post-hoc tests
  "performance",  # Diagnostics and assumptions
  "vegan",        # NMDS, diversity
  "data.table",   # Fast data manipulation
  "ggpubr",       # Publication plots
  "patchwork",    # Plot combining
  "effectsize",   # eta² etc
  "gptstudio",    # RStudio ChatGPT Addin
  "reticulate",   # Python + AI packages in R
  "chatgpt"       # OpenAI API wrapper (if using GPT via R)
))

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(vegan)
library(data.table)
library(ggpubr)
library(patchwork)
library(effectsize)

# Optional AI tools
library(gptstudio)   # AI integration in RStudio
library(chatgpt)     # OpenAI wrapper (needs API key)
library(reticulate)  # Python + R combo

install.packages("emmeans")




install.packages("tidyverse")

setwd("C:/Users/agube/OneDrive/Documentos/r_programming")

library(tidyverse)
library(readxl)

# Load 5 raw field files
processed     <- read_excel("2025_field_processed_samples.xlsx")
initial_mass  <- read_excel("2025_field_initial_mass.xlsx")
transport     <- read_excel("2025_field_transport_controls.xlsx")
site_info     <- read_excel("2025_field_site_information_Kinzig.xlsx")
inverts       <- read_excel("2025_field_invertebrates_coarse_leaf_bags.xlsx")

site_info <- site_info %>%
  mutate(
    LandUse = case_when(
      tolower(trimws(rural)) == "x" ~ "rural",
      tolower(trimws(urban)) == "x" ~ "urban",
      TRUE ~ NA_character_
    )
  ) %>%
  select(stream, LandUse)


# Calculate mean transport control values per stream
transport_summary <- transport %>%
  group_by(stream) %>%
  summarise(
    mean_OMD60 = mean(OMD_60_mg, na.rm = TRUE),
    mean_OMD500 = mean(OMD_500_mg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    AFDM_0_corrected = mean_OMD60 - mean_OMD500
  )
colnames(transport)

transport <- transport %>%
  mutate(AFDM_0 = mass_60_g - mass_500_g)

mean_AFDM_0 <- mean(transport$AFDM_0, na.rm = TRUE)


# Filter Kinzig catchment
kinzig_processed <- processed %>%
  filter(catchment == "Kinzig") %>%
  mutate(
    AFDM_final = OMD_60_mg - OMD_500_mg,
    AFDM_0 = mean_AFDM_0,  # apply same corrected AFDM₀ to all
    deg_days = 14 * 12.5,  # change if needed
    OMD_dday = ((AFDM_0 - AFDM_final) / AFDM_0) * 100 / deg_days
  )

# Filter Kinzig catchment
kinzig_processed <- processed %>%
  filter(catchment == "Kinzig") %>%
  mutate(
    AFDM_final = OMD_60_mg - OMD_500_mg,
    AFDM_0 = mean_AFDM_0,  # apply same corrected AFDM₀ to all
    deg_days = 14 * 12.5,  # change if needed
    OMD_dday = ((AFDM_0 - AFDM_final) / AFDM_0) * 100 / deg_days
  )

processed <- processed %>%
  mutate(
    OMD_60_mg = as.numeric(OMD_60_mg),
    OMD_500_mg = as.numeric(OMD_500_mg)
  )

processed %>%
  summarise(
    missing_60 = sum(is.na(OMD_60_mg)),
    missing_500 = sum(is.na(OMD_500_mg))
  )

# Join land use info
kinzig_labeled <- kinzig_processed %>%
  left_join(site_info, by = "stream")

# Check and clean
kinzig_labeled <- kinzig_labeled %>%
  filter(!is.na(LandUse), !is.na(OMD_dday))  # keep only usable rows

# Reshape wide: one row per replicate, with fine & coarse OMD side-by-side
model_data <- kinzig_labeled %>%
  select(stream, replicate, type, LandUse, OMD_dday) %>%
  pivot_wider(names_from = type, values_from = OMD_dday) %>%
  rename(
    k_micro = fine,
    k_total = coarse
  ) %>%
  mutate(
    k_macro = k_total - k_micro
  )

processed <- processed %>%
  mutate(
    OMD_60_mg = as.numeric(OMD_60_mg),
    OMD_500_mg = as.numeric(OMD_500_mg)
  )

# Rebuild kinzig_processed safely
kinzig_processed <- processed %>%
  filter(catchment == "Kinzig") %>%
  mutate(
    AFDM_final = OMD_60_mg - OMD_500_mg,
    AFDM_0 = mean_AFDM_0,
    deg_days = 14 * 12.5,
    OMD_dday = ((AFDM_0 - AFDM_final) / AFDM_0) * 100 / deg_days
  )



# Add land use
kinzig_labeled <- kinzig_processed %>%
  left_join(site_info, by = "stream") %>%
  filter(!is.na(LandUse), !is.na(OMD_dday))

# Pivot to wide format
model_data <- kinzig_labeled %>%
  select(stream, replicate, type, LandUse, OMD_dday) %>%
  pivot_wider(names_from = type, values_from = OMD_dday) %>%
  rename(
    k_micro = fine,
    k_total = coarse
  ) %>%
  mutate(
    k_macro = k_total - k_micro
  )

#
  

glimpse(model_data)


library(lme4)
library(lmerTest)
library(emmeans)

model_micro <- lmer(k_micro ~ LandUse + (1 | stream), data = model_data)
summary(model_micro)
emmeans(model_micro, pairwise ~ LandUse)


model_total <- lmer(k_total ~ LandUse + (1 | stream), data = model_data)
summary(model_total)
emmeans(model_total, pairwise ~ LandUse)


model_macro <- lmer(k_macro ~ LandUse + (1 | stream), data = model_data)
summary(model_macro)
emmeans(model_macro, pairwise ~ LandUse)


model_micro <- lmer(k_micro ~ LandUse + (1 | stream), data = model_data)

# Summary (estimates, p-values)
summary(model_micro)

# Post-hoc pairwise contrasts (emmeans)
emmeans(model_micro, pairwise ~ LandUse)


lmer(k_micro ~ LandUse + (1 | stream))

colnames(model_data)

library(lme4)
library(lmerTest)
library(emmeans)

model_micro <- lmer(k_micro ~ LandUse + (1 | stream), data = model_data)

# Summary
summary(model_micro)

# Post-hoc contrast
emmeans(model_micro, pairwise ~ LandUse)


# ------------------------------------------------------------
# Hypothesis 2: Total decomposition is higher in rural streams
# Test whether leaf litter decomposition from coarse mesh bags 
# (including microbial + macroinvertebrate contributions) differs 
# between rural and urban streams in the Kinzig catchment.
# Linear mixed-effects model with stream as random intercept.
# ------------------------------------------------------------


model_total <- lmer(k_total ~ LandUse + (1 | stream), data = model_data)
summary(model_total)
emmeans(model_total, pairwise ~ LandUse)


# ------------------------------------------------------------
# Hypothesis 2: Total decomposition is higher in rural streams
# Model: k_total ~ LandUse + (1 | stream)
# 
# Results:
# Estimated marginal means:
#   rural: -0.81 ± 0.34
#   urban: -1.56 ± 0.35
# Contrast (rural - urban): +0.75, p = 0.142 (not significant)
# 
# Interpretation:
# Total decomposition was higher in rural streams, but the 
# difference was not statistically significant. The estimated 
# effect size suggests a potential trend, but with considerable 
# variability across streams.
# ------------------------------------------------------------



# ------------------------------------------------------------
# Hypothesis 3: Macroinvertebrate contribution is higher in rural streams
# Model: k_macro ~ LandUse + (1 | stream)
#
# Test whether the difference between coarse and fine mesh bags 
# (representing invertebrate-mediated decomposition) is greater 
# in rural than in urban streams in the Kinzig catchment.
# ------------------------------------------------------------
 

k_macro ~ LandUse + (1 | stream)


model_macro <- lmer(k_macro ~ LandUse + (1 | stream), data = model_data)
summary(model_macro)
emmeans(model_macro, pairwise ~ LandUse)


# ------------------------------------------------------------
# Hypothesis 3: Macroinvertebrate contribution is higher in rural streams
# Model: k_macro ~ LandUse + (1 | stream)
#
# Results:
# Estimated marginal means:
#   rural: 1.66 ± 0.34
#   urban: 0.87 ± 0.36
# Contrast (rural - urban): +0.79, p = 0.128 (not significant)
#
# Interpretation:
# Macroinvertebrate-mediated decomposition tended to be higher in 
# rural streams, aligning with the hypothesis. However, the difference 
# was not statistically significant. The effect size suggests a trend, 
# but inter-stream variation remains high.
# ------------------------------------------------------------


colnames(inverts)


# Step 1: Keep only coarse mesh invertebrate data
inverts_clean <- inverts %>%
  filter(type == "coarse") %>%
  mutate(
    total_abundance = rowSums(across(c(Gammarids, Trichoptera, Chironomids, Others)), na.rm = TRUE)
  ) %>%
  select(stream, replicate, total_abundance, Gammarids, Trichoptera, Chironomids, Others)

# Step 2: Join with your model_data
model_data <- model_data %>%
  left_join(inverts_clean, by = c("stream", "replicate"))


# Step 1: Filter to coarse bags and convert counts to numeric
inverts_clean <- inverts %>%
  filter(type == "coarse") %>%
  mutate(
    Gammarids    = as.numeric(Gammarids),
    Trichoptera  = as.numeric(Trichoptera),
    Chironomids  = as.numeric(Chironomids),
    Others       = as.numeric(Others),
    total_abundance = rowSums(across(c(Gammarids, Trichoptera, Chironomids, Others)), na.rm = TRUE)
  ) %>%
  select(stream, replicate, total_abundance, Gammarids, Trichoptera, Chironomids, Others)


glimpse(model_data)


# Fix spacing in stream and replicate columns
inverts_clean <- inverts_clean %>%
  mutate(
    stream = str_trim(stream),
    replicate = str_trim(replicate)
  )

model_data <- model_data %>%
  mutate(
    stream = str_trim(stream),
    replicate = str_trim(replicate)
  )

# Re-run the join
model_data <- model_data %>%
  left_join(inverts_clean, by = c("stream", "replicate"))


glimpse(model_data)

ggplot(model_data, aes(x = total_abundance, y = k_macro)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(x = "Total Invertebrate Abundance", y = "Macroinvertebrate-Mediated Decomposition (k_macro)",
       title = "Relationship Between Invertebrate Abundance and Leaf Breakdown")



ggplot(model_data, aes(x = total_abundance, y = k_macro)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_bw() +
  labs(
    title = "Relationship Between Invertebrate Abundance and k_macro",
    x = "Total Macroinvertebrate Abundance",
    y = "Macroinvertebrate-Mediated Decomposition (k_macro)"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )


model_macro_abund <- lmer(k_macro ~ LandUse + total_abundance + (1 | stream), data = model_data)

# Summary with estimates
summary(model_macro_abund)

# Post-hoc for LandUse
emmeans(model_macro_abund, pairwise ~ LandUse)


# ------------------------------------------------------------
# Figure 1: Boxplots of decomposition rates by Land Use
# This figure compares k_micro, k_total, and k_macro between 
# urban and rural streams using boxplots. It visualizes raw 
# distribution across sites without model smoothing.
# ------------------------------------------------------------

# Prepare tidy long format for plotting
plot_data <- model_data %>%
  pivot_longer(cols = c(k_micro, k_total, k_macro),
               names_to = "Decomposition_Type",
               values_to = "k_value")

# Create the plot
ggplot(plot_data, aes(x = LandUse, y = k_value, fill = LandUse)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  facet_wrap(~ Decomposition_Type, scales = "free_y") +
  labs(
    title = "Figure 1. Decomposition Rates by Land Use",
    y = "Decomposition Rate (per degree day)",
    x = "Land Use"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )


# Color palette: consistent with Schreiner (e.g. rural = green, urban = orange)
ggplot(plot_data, aes(x = LandUse, y = k_value, fill = LandUse)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  facet_wrap(~ Decomposition_Type, scales = "free_y") +
  scale_fill_manual(
    values = c("rural" = "forestgreen", "urban" = "darkorange"),
    labels = c("Rural", "Urban")
  ) +
  labs(
    title = "Figure 1. Decomposition Rates by Land Use",
    y = "Decomposition Rate (per degree day)",
    x = "Land Use",
    fill = "Land Use"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# ------------------------------------------------------------
# Figure 2: Estimated Marginal Means from Mixed Models
# Displays modeled mean decomposition rates for rural and 
# urban streams across three response variables (k_micro, 
# k_total, k_macro). Error bars show 95% confidence intervals.
# ------------------------------------------------------------


# Get emmeans for all 3 response variables
emm_micro <- emmeans(model_micro, ~ LandUse) %>% as.data.frame() %>% mutate(type = "k_micro")
emm_total <- emmeans(model_total, ~ LandUse) %>% as.data.frame() %>% mutate(type = "k_total")
emm_macro <- emmeans(model_macro, ~ LandUse) %>% as.data.frame() %>% mutate(type = "k_macro")

# Combine all
emm_all <- bind_rows(emm_micro, emm_total, emm_macro)

# Plot
ggplot(emm_all, aes(x = LandUse, y = emmean, fill = LandUse)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_wrap(~ type, scales = "free_y") +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(
    title = "Figure 2. Estimated Marginal Means of Decomposition Rates",
    x = "Land Use",
    y = "Estimated k (per degree day)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )


# ------------------------------------------------------------
# Figure 3: Relationship Between Invertebrate Abundance and k_macro
# Scatterplot showing how total macroinvertebrate abundance relates 
# to macroinvertebrate-mediated decomposition. Includes linear model 
# fit with 95% confidence interval and points colored by Land Use.
# ------------------------------------------------------------

ggplot(model_data, aes(x = total_abundance, y = k_macro, color = LandUse)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(
    title = "Figure 3. Invertebrate Abundance vs Macroinvertebrate-Mediated Decomposition",
    x = "Total Macroinvertebrate Abundance",
    y = "k_macro (coarse − fine decomposition)",
    color = "Land Use"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )



# ------------------------------------------------------------
# Table 1: Summary of Fixed Effects from Mixed Models
# Extract fixed effect estimates, standard errors, degrees of freedom, 
# t-values and p-values for each model and combine into a table.
# ------------------------------------------------------------


# Extract fixed effects
get_fixed <- function(model, response_name) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    mutate(response = response_name) %>%
    select(response, term, estimate, std.error, df, statistic, p.value)
}

# Load broom.mixed (if not yet installed)
if (!requireNamespace("broom.mixed", quietly = TRUE)) install.packages("broom.mixed")
library(broom.mixed)

# Get rows from each model
fx_micro <- get_fixed(model_micro, "k_micro")
fx_total <- get_fixed(model_total, "k_total")
fx_macro <- get_fixed(model_macro, "k_macro")
fx_macro_abund <- get_fixed(model_macro_abund, "k_macro (+abundance)")

# Combine
table1 <- bind_rows(fx_micro, fx_total, fx_macro, fx_macro_abund)

# Round for reporting
table1_clean <- table1 %>%
  mutate(across(c(estimate, std.error, statistic, p.value), ~ round(.x, 3)))

# Show preview
print(table1_clean)

# ------------------------------------------------------------
# Table 3: Pairwise Contrasts Between Rural and Urban Streams
# Extracted from emmeans pairwise comparisons for each response
# variable. Includes estimate, SE, df, t ratio, and p-value.
# ------------------------------------------------------------


# Get pairwise contrasts for each model
cont_micro <- emmeans(model_micro, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_micro")

cont_total <- emmeans(model_total, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_total")

cont_macro <- emmeans(model_macro, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_macro")

cont_macro_abund <- emmeans(model_macro_abund, pairwise ~ LandUse)$contrasts %>%
  as.data.frame() %>% mutate(response = "k_macro (+abundance)")

# Combine
table3 <- bind_rows(cont_micro, cont_total, cont_macro, cont_macro_abund) %>%
  select(response, contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(across(c(estimate, SE, t.ratio, p.value), ~ round(.x, 3)))

# Show result
print(table3)


# ------------------------------------------------------------
# Figure 4: Random Stream Intercepts from Mixed Model
# Shows stream-level deviations from the global intercept for 
# k_macro model. Visualizes variability captured by the random 
# effect term (1 | stream).
# ------------------------------------------------------------


# Extract stream-level intercepts from model
stream_intercepts <- ranef(model_macro)$stream %>%
  as.data.frame() %>%
  rename(deviation = `(Intercept)`) %>%
  rownames_to_column(var = "stream")

# Plot
ggplot(stream_intercepts, aes(x = reorder(stream, deviation), y = deviation)) +
  geom_point(size = 3, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Figure 4. Stream-wise Random Intercepts for k_macro",
    y = "Deviation from Overall Mean",
    x = "Stream"
  )


# ------------------------------------------------------------
# Supplemental Calculations: Effect Sizes and Rate Conversions
# Calculates Cohen's d, percent difference (rural vs urban),
# and decomposition rate per day (instead of per degree day).
# Uses emmeans outputs for k_micro, k_total, and k_macro.
# ------------------------------------------------------------


# Rebuild emmeans table manually (or reuse from Table 2)
emm_summary <- tibble::tibble(
  response = c("k_micro", "k_total", "k_macro"),
  rural_mean = c(-2.47, -0.811, 1.656),
  rural_se = c(0.086, 0.337, 0.342),
  urban_mean = c(-2.43, -1.562, 0.865),
  urban_se = c(0.089, 0.352, 0.357),
  temp_C = 12.5
)

# Calculate:
# - Cohen's d
# - Percent difference
# - Per day decomposition rates
calc_metrics <- emm_summary %>%
  mutate(
    cohen_d = (rural_mean - urban_mean) / ((rural_se + urban_se) / 2),
    percent_diff = ((rural_mean - urban_mean) / abs(urban_mean)) * 100,
    rural_per_day = rural_mean * temp_C,
    urban_per_day = urban_mean * temp_C
  )

# Show results
print(calc_metrics)

write.csv(calc_metrics, "table4_supplemental_calculations.csv", row.names = FALSE)


# ------------------------------------------------------------
# Annotated Figure 2: Estimated Marginal Means with % Differences
# Adds numeric annotations showing percent difference between 
# rural and urban means for each decomposition type.
# ------------------------------------------------------------

# Add annotation manually (from Table 4)
annotations <- tibble::tibble(
  type = c("k_micro", "k_total", "k_macro"),
  label = c("-1.6%", "+48.1%", "+91.4%"),  # from percent_diff
  y_pos = c(-2.2, -0.2, 2.3)  # adjust above the taller bar
)

# Plot again with annotation layer
ggplot(emm_all, aes(x = LandUse, y = emmean, fill = LandUse)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_wrap(~ type, scales = "free_y") +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  geom_text(data = annotations, aes(x = 1.5, y = y_pos, label = label), inherit.aes = FALSE,
            size = 4.5, fontface = "bold") +
  labs(
    title = "Figure 2. Estimated Marginal Means of Decomposition Rates",
    x = "Land Use", y = "Estimated k (per degree day)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

library(performance)

# Diagnostics for each core model
check_model(model_micro)
check_model(model_total)
check_model(model_macro)
check_model(model_macro_abund)  # optional

# Optional: Save residuals
resids <- resid(model_macro)
fitted_vals <- fitted(model_macro)

# Simple residual plot
plot(fitted_vals, resids,
     main = "Residuals vs Fitted (k_macro)",
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)


# ------------------------------------------------------------
# Step 1: Prepare Physico-Chemical Data for Merging
# Clean column names and reduce to one row per stream.
# ------------------------------------------------------------

library(readxl)
library(janitor)
library(dplyr)

# Load the Excel file
phch_raw <- read_excel("202505_phch_data_Kinzig.xlsx", sheet = 1) %>%
  clean_names()

# Preview
glimpse(phch_raw)

# Select and rename key columns
phch_clean <- phch_raw %>%
  select(stream,
         ph = p_h,
         water_temp = water_temp,
         oxygen_dissolved = oxygen_dissolved,
         oxygen_saturation = oxygen_saturation,
         redox = redox_potential,
         conductivity = conduc_tivity,
         water_level = water_level_at_pole) %>%
  filter(!is.na(stream)) %>%
  mutate(stream = str_trim(stream))


# ------------------------------------------------------------
# Step 2: Merge Physico-Chemical Data with model_data
# Join by stream and retain all decomposition + environmental variables.
# ------------------------------------------------------------

# Trim stream names in model_data just in case
model_data <- model_data %>%
  mutate(stream = str_trim(stream))

# Join
model_data <- model_data %>%
  left_join(phch_clean, by = "stream")

# Check result
glimpse(model_data)


k_macro ~ LandUse + water_temp + (1 | stream)

# ------------------------------------------------------------
# Step 3: Model Effect of Water Temperature on k_macro
# Tests whether stream temperature explains decomposition
# differences beyond land use effects.
# ------------------------------------------------------------

model_temp <- lmer(k_macro ~ LandUse + water_temp + (1 | stream), data = model_data)

summary(model_temp)
emmeans(model_temp, pairwise ~ LandUse)



# ------------------------------------------------------------
# Figure 5: Relationship Between Water Temperature and k_macro
# Scatterplot showing macroinvertebrate-mediated decomposition
# as a function of mean stream temperature. Colored by Land Use.
# ------------------------------------------------------------


ggplot(model_data, aes(x = water_temp, y = k_macro, color = LandUse)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  labs(
    title = "Figure 5. Water Temperature vs k_macro",
    x = "Water Temperature (°C)",
    y = "k_macro (coarse − fine decomposition)",
    color = "Land Use"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )



# ------------------------------------------------------------
# Combined Environmental Model: Water Temp, Conductivity, Oxygen
# Tests the additive effect of stream physico-chemical variables
# on macroinvertebrate-mediated decomposition.
# ------------------------------------------------------------

model_env <- lmer(k_macro ~ LandUse + water_temp + conductivity + oxygen_dissolved + (1 | stream),
                  data = model_data)

summary(model_env)



# ------------------------------------------------------------
# Figure 6: Partial Effects of Environmental Variables on k_macro
# Shows predicted relationships from the additive model for 
# temperature, conductivity, and dissolved oxygen.
# Adjusted for other predictors (partial effects).
# ------------------------------------------------------------


library(ggeffects)
library(patchwork)

# Get partial predictions
eff_temp <- ggpredict(model_env, terms = "water_temp")
eff_cond <- ggpredict(model_env, terms = "conductivity")
eff_oxy  <- ggpredict(model_env, terms = "oxygen_dissolved")

# Build individual plots
p1 <- plot(eff_temp) + labs(title = "Water Temperature (°C)", y = "Predicted k_macro") +
  theme_bw()

p2 <- plot(eff_cond) + labs(title = "Conductivity (µS/cm)", y = "") +
  theme_bw()

p3 <- plot(eff_oxy) + labs(title = "Dissolved Oxygen (mg/L)", y = "") +
  theme_bw()

# Combine
(p1 | p2 | p3) +
  plot_annotation(title = "Figure 6. Partial Effects of Environmental Predictors on Macroinvertebrate Decomposition")


install.packages("patchwork")
installed.packages(ggeffects)
install.packages("ggeffects")
library(ggeffects)

eff_temp <- ggpredict(model_env, terms = "water_temp")
eff_cond <- ggpredict(model_env, terms = "conductivity")
eff_oxy  <- ggpredict(model_env, terms = "oxygen_dissolved")

library(patchwork)

p1 <- plot(eff_temp) + labs(title = "Water Temperature (°C)", y = "Predicted k_macro") + theme_bw()
p2 <- plot(eff_cond) + labs(title = "Conductivity (µS/cm)", y = "") + theme_bw()
p3 <- plot(eff_oxy) + labs(title = "Dissolved Oxygen (mg/L)", y = "") + theme_bw()

(p1 | p2 | p3) +
  plot_annotation(title = "Figure 6. Partial Effects of Environmental Predictors on Macroinvertebrate Decomposition")


library(ggplot2)


p1 <- plot(eff_temp) + labs(title = "Water Temperature (°C)", y = "Predicted k_macro") + theme_bw()
p2 <- plot(eff_cond) + labs(title = "Conductivity (µS/cm)", y = "") + theme_bw()
p3 <- plot(eff_oxy) + labs(title = "Dissolved Oxygen (mg/L)", y = "") + theme_bw()

(p1 | p2 | p3) +
  plot_annotation(title = "Figure 6. Partial Effects of Environmental Predictors on Macroinvertebrate Decomposition")


