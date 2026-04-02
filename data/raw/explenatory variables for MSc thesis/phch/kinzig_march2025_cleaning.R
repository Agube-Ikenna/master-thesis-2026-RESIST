# =============================================================================
# Kinzig River – March 2025 Data Cleaning
# =============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(writexl)

options(scipen = 999)  # disable scientific notation

kinzig_path <- "SFB_FieldStudies_PhysicoChemicalData_Kinzig_2026-01-28.xlsx"

# -----------------------------------------------------------------------------
# Helper: read a sheet safely (handles mismatched header/data column counts)
# -----------------------------------------------------------------------------
read_kinzig_sheet <- function(path, sheet_name) {
  df_raw        <- read_excel(path, sheet = sheet_name, skip = 3, col_names = FALSE)
  n_cols        <- ncol(df_raw)
  header_row    <- read_excel(path, sheet = sheet_name, col_names = FALSE, n_max = 1)
  col_names_vec <- as.character(header_row[1, 1:n_cols])
  col_names_vec <- ifelse(
    is.na(col_names_vec) | col_names_vec == "" | col_names_vec == "NA",
    paste0("col_", seq_along(col_names_vec)),
    col_names_vec
  )
  col_names_vec    <- make.names(col_names_vec, unique = TRUE)
  names(df_raw)    <- col_names_vec
  names(df_raw)[1] <- "date"
  df_raw %>%
    mutate(date = as.Date(date)) %>%
    filter(!is.na(date), year(date) == 2025, month(date) == 3)
}

# =============================================================================
# 1. READ & FILTER TO MARCH 2025
# =============================================================================
kinzig_field_mar25        <- read_kinzig_sheet(kinzig_path, "field | bi-weekly")
kinzig_lab_biweekly_mar25 <- read_kinzig_sheet(kinzig_path, "lab | bi-weekly")
kinzig_lab_monthly_mar25  <- read_kinzig_sheet(kinzig_path, "lab | once a month")

cat("Rows after filtering to March 2025:\n")
cat("  field | bi-weekly :", nrow(kinzig_field_mar25), "\n")
cat("  lab   | bi-weekly :", nrow(kinzig_lab_biweekly_mar25), "\n")
cat("  lab   | once/month:", nrow(kinzig_lab_monthly_mar25), "\n\n")

# =============================================================================
# 2. AVERAGE BY SITE + RENAME COLUMNS + ROUND DECIMALS (all in one step)
# =============================================================================

# Detect exact column names (may vary in dot count after make.names)
ph_col   <- grep("^pH",          names(kinzig_field_mar25), value = TRUE)[1]
temp_col <- grep("water.temp",   names(kinzig_field_mar25), value = TRUE)[1]

cat("Detected pH column       :", ph_col,   "\n")
cat("Detected water temp column:", temp_col, "\n")

# ── field | bi-weekly ─────────────────────────────────────────────────────────
kinzig_field_merged <- kinzig_field_mar25 %>%
  group_by(site = site.name) %>%
  summarise(
    pH_dimensionless      = round(mean(as.numeric(.data[[ph_col]]),   na.rm = TRUE), 2),
    conductivity_uS_cm    = round(mean(conduc.tivity,                 na.rm = TRUE), 1),
    dissolved_oxygen_mg_l = round(mean(dissolved.oxygen,              na.rm = TRUE), 2),
    water_temp_degC       = round(mean(as.numeric(.data[[temp_col]]), na.rm = TRUE), 2),
    water_level_cm        = round(mean(water.level.at.pole,           na.rm = TRUE), 1),
    .groups = "drop"
  )

# ── lab | bi-weekly ───────────────────────────────────────────────────────────
kinzig_lab_biweekly_merged <- kinzig_lab_biweekly_mar25 %>%
  group_by(site = site.name) %>%
  summarise(
    ortho_phosphate_mg_l = round(mean(ortho.phosphate, na.rm = TRUE), 4),
    total_phosphate_mg_l = round(mean(total.phosphate, na.rm = TRUE), 4),
    ammonia_mg_l         = round(mean(ammonia,         na.rm = TRUE), 4),
    nitrite_mg_l         = round(mean(nitrite,         na.rm = TRUE), 4),
    nitrate_mg_l         = round(mean(nitrate,         na.rm = TRUE), 3),
    total_nitrogen_mg_l  = round(mean(total.nitrogen,  na.rm = TRUE), 4),
    chloride_mg_l        = round(mean(chloride,        na.rm = TRUE), 3),
    sulfate_mg_l         = round(mean(sulfate,         na.rm = TRUE), 3),
    .groups = "drop"
  )

# ── lab | once a month ────────────────────────────────────────────────────────
kinzig_lab_monthly_merged <- kinzig_lab_monthly_mar25 %>%
  group_by(site = site.name) %>%
  summarise(
    Cd_mg_l   = round(mean(Cd,                na.rm = TRUE), 6),
    Cu_mg_l   = round(mean(Cu,                na.rm = TRUE), 6),
    Fe_mg_l   = round(mean(Fe,                na.rm = TRUE), 4),
    Ni_mg_l   = round(mean(Ni,                na.rm = TRUE), 6),
    Pb_mg_l   = round(mean(Pb,                na.rm = TRUE), 6),
    Zn_mg_l   = round(mean(Zn,                na.rm = TRUE), 6),
    Ca_mg_l   = round(mean(Ca2.,              na.rm = TRUE), 3),
    Mg_mg_l   = round(mean(Mg2.,              na.rm = TRUE), 3),
    Na_mg_l   = round(mean(Na.,               na.rm = TRUE), 3),
    K_mg_l    = round(mean(K.,                na.rm = TRUE), 3),
    Co3_mg_l  = round(mean(Co32.,             na.rm = TRUE), 6),
    HCO3_mg_l = round(mean(as.numeric(HCO3.), na.rm = TRUE), 3),
    .groups = "drop"
  )

# =============================================================================
# 3. JOIN ALL THREE SHEETS INTO ONE (by site)
# =============================================================================
kinzig_mar2025 <- kinzig_field_merged %>%
  full_join(kinzig_lab_biweekly_merged, by = "site") %>%
  full_join(kinzig_lab_monthly_merged,  by = "site")

# =============================================================================
# 4. FINAL CHECK
# =============================================================================
cat("Final joined rows:", nrow(kinzig_mar2025), "\n")
cat("Final joined cols:", ncol(kinzig_mar2025), "\n")
print(names(kinzig_mar2025))
print(kinzig_mar2025)

# =============================================================================
# 5. SAVE TO EXCEL
# =============================================================================
write_xlsx(kinzig_mar2025, "kinzig_mar2025_final.xlsx")
cat("\nSaved: kinzig_mar2025_final.xlsx\n")
