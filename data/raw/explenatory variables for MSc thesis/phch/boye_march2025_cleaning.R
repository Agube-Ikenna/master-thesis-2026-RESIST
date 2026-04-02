# =============================================================================
# Boye River – March 2025 Data Cleaning
# =============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(writexl)

options(scipen = 999)  # disable scientific notation

boye_path <- "SFB_FieldStudies_PhysicoChemicalData_Boye_2026-01-28.xlsx"

# -----------------------------------------------------------------------------
# Helper: read a sheet safely (handles mismatched header/data column counts)
# -----------------------------------------------------------------------------
read_boye_sheet <- function(path, sheet_name) {
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
boye_field_mar25        <- read_boye_sheet(boye_path, "field | bi-weekly")
boye_lab_biweekly_mar25 <- read_boye_sheet(boye_path, "lab | bi-weekly")
boye_lab_monthly_mar25  <- read_boye_sheet(boye_path, "lab | once a month")

cat("Rows after filtering to March 2025:\n")
cat("  field | bi-weekly :", nrow(boye_field_mar25), "\n")
cat("  lab   | bi-weekly :", nrow(boye_lab_biweekly_mar25), "\n")
cat("  lab   | once/month:", nrow(boye_lab_monthly_mar25), "\n\n")

# =============================================================================
# 2. AVERAGE BY SITE + RENAME COLUMNS + ROUND DECIMALS (all in one step)
# =============================================================================

# Detect exact oxygen saturation column name (may vary: one or two dots)
oxy_col <- grep("oxygen.*saturation", names(boye_field_mar25), value = TRUE)[1]
cat("Detected oxygen saturation column:", oxy_col, "\n")

# ── field | bi-weekly ─────────────────────────────────────────────────────────
boye_field_merged <- boye_field_mar25 %>%
  group_by(site = site.name) %>%
  summarise(
    pH_dimensionless      = round(mean(pH,                                  na.rm = TRUE), 2),
    conductivity_uS_cm    = round(mean(conduc.tivity,                       na.rm = TRUE), 1),
    dissolved_oxygen_mg_l = round(mean(dissolved.oxygen,                    na.rm = TRUE), 2),
    oxygen_saturation_pct = round(mean(as.numeric(.data[[oxy_col]]),         na.rm = TRUE), 1),
    water_temp_degC       = round(mean(water.temp,                          na.rm = TRUE), 2),
    water_level_cm        = round(mean(water.level.at.pole,                 na.rm = TRUE), 1),
    .groups = "drop"
  )

# ── lab | bi-weekly ───────────────────────────────────────────────────────────
boye_lab_biweekly_merged <- boye_lab_biweekly_mar25 %>%
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
boye_lab_monthly_merged <- boye_lab_monthly_mar25 %>%
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
# 3. FINAL CHECK
# =============================================================================
cat("Final column names:\n")
cat("  field | bi-weekly :\n");  print(names(boye_field_merged))
cat("  lab   | bi-weekly :\n");  print(names(boye_lab_biweekly_merged))
cat("  lab   | once/month:\n");  print(names(boye_lab_monthly_merged))

cat("\nFinal row counts:\n")
cat("  field | bi-weekly :", nrow(boye_field_merged), "\n")
cat("  lab   | bi-weekly :", nrow(boye_lab_biweekly_merged), "\n")
cat("  lab   | once/month:", nrow(boye_lab_monthly_merged), "\n")

cat("\nPreview:\n")
print(boye_field_merged)
print(boye_lab_biweekly_merged)
print(boye_lab_monthly_merged)

# =============================================================================
# 4. SAVE TO EXCEL
# =============================================================================
write_xlsx(
  list(
    "field_biweekly" = boye_field_merged,
    "lab_biweekly"   = boye_lab_biweekly_merged,
    "lab_monthly"    = boye_lab_monthly_merged
  ),
  "boye_mar2025_cleaned.xlsx"
)

cat("\nSaved: boye_mar2025_cleaned.xlsx\n")
