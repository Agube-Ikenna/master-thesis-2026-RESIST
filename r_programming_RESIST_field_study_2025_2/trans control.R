# --- 1. Load Data ---
transport <- read.csv("2025_field_transport_controls.csv", na.strings = c("NR", "NA"))
initial_mass <- read.csv("2025_field_initial_mass.csv", na.strings = c("NR", "NA"))
processed_kinzig <- read.csv("2025_field_processed_samples.csv", na.strings = c("NR", "NA"))

# --- 2. Standardize Key Columns ---
transport$type <- tolower(trimws(transport$type))
transport$catchment <- tolower(trimws(transport$catchment))
initial_mass$type <- tolower(trimws(initial_mass$type))
processed_kinzig$type <- tolower(trimws(processed_kinzig$type))

# --- 3. Calculate Handling Loss Ratio and Flag Outliers ---
transport$perc_60 <- transport$mass_60_g / transport$initial_mass_g
# You can set a more refined rule or use supervisor notes:
transport$outlier <- with(transport, perc_60 < 0.8 | perc_60 > 1.05)  # adjust as needed

# --- 4. Filter to Kinzig and Remove Outliers ---
transport_kinzig <- transport %>%
  filter(catchment == "kinzig" & !outlier)

# --- 5. Calculate Correction Factors by Mesh Type ---
correction_factors_kinzig <- transport_kinzig %>%
  group_by(type) %>%
  summarise(perc_60 = mean(perc_60, na.rm = TRUE))

print(correction_factors_kinzig)
# # A tibble: 2 × 2
#   type   perc_60
#   <chr>    <dbl>
# 1 coarse   0.897
# 2 fine     0.834

# --- 6. Merge Initial Mass if Needed ---
if(!"initial_mass_g" %in% names(processed_kinzig)) {
  processed_kinzig <- processed_kinzig %>%
    left_join(initial_mass[, c("bag_ID", "type", "initial_mass_g")], by = c("bag_ID", "type"))
}

# --- 7. Join Correction Factors and Calculate AFDM_start ---
processed_kinzig <- processed_kinzig %>%
  left_join(correction_factors_kinzig, by = "type") %>%
  mutate(
    AFDM_start = initial_mass_g * perc_60
  )

# --- 8. QA: Check Output ---
summary(processed_kinzig$AFDM_start)
head(processed_kinzig[, c("bag_ID", "type", "initial_mass_g", "perc_60", "AFDM_start")])
