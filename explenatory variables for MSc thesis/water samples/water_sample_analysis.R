library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

# ── File paths ────────────────────────────────────────────────────────────────
boye_path   <- "2025Wasserprobenanalyse_SFB_RESIST_Boye_System.xlsx"
kinzig_path <- "2025Wasserprobenanalyse_SFB_RESIST_Kinzig_System.xlsx"


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

drop_junk <- function(df) {
  df |> select(-matches("^\\.\\.\\.\\d+"))
}

drop_empty_cols <- function(df) {
  df |> select(where(~ !all(is.na(.))))
}

force_numeric <- function(df, id_cols = c("date", "site", "system")) {
  df |> mutate(across(-any_of(id_cols), ~ round(as.numeric(.), 4)))
}

to_long <- function(df, id_cols = c("date", "site", "system")) {
  df |>
    pivot_longer(cols      = -all_of(id_cols),
                 names_to  = "parameter",
                 values_to = "value") |>
    filter(!is.na(value))
}


# ══════════════════════════════════════════════════════════════════════════════
# 1.  READ, FILTER TO MARCH 2025, CLEAN DATE
# ══════════════════════════════════════════════════════════════════════════════

read_march <- function(path, sheet, date_col) {
  read_excel(path, sheet = sheet) |>
    rename(date = all_of(date_col)) |>
    filter(format(date, "%Y-%m") == "2025-03") |>
    mutate(date = as.Date(date)) |>        # strip time, keep date only
    drop_junk() |>
    drop_empty_cols()
}

# ── Boye ----------------------------------------------------------------------
boye_freq <- read_march(boye_path, "every watersample", "Datum") |>
  rename(
    site            = Probestelle,
    ortho_phosphate = `ortho-phosphate`,
    total_phosphate = `total phosphate mg/l`,
    ammonia         = ammonia,
    nitrite         = `nitrite mg/l`,
    nitrate         = `nitrate mg/l`,
    total_nitrogen  = `total nitrogen`,
    chloride        = `chloride mg/l`,
    sulfate         = `sulfate mg/l`
  ) |>
  select(-any_of("CODE")) |>
  mutate(system = "Boye") |>
  force_numeric()

boye_month <- read_march(boye_path, "once a month", "Datum") |>
  rename(
    site      = Probestelle,
    Cd        = `Cd mg/l`,
    Cu        = `Cu mg/l`,
    Fe        = `Fe mg/l`,
    Ni        = `Ni mg/l`,
    Pb        = `Pb mg/l`,
    Zn        = `Zn mg/l`,
    Ca        = `Ca2+`,
    Mg        = `Mg2+`,
    Na        = `Na+`,
    K         = `K+`,
    carbonate = `Co32-`,
    HCO3      = `HCO3- mg/l`
  ) |>
  select(-any_of("CODE")) |>
  mutate(system = "Boye") |>
  force_numeric()

# ── Kinzig --------------------------------------------------------------------
kinzig_freq <- read_march(kinzig_path, "every watersample", "date") |>
  rename(
    ortho_phosphate = `ortho-phosphate mg/l`,
    total_phosphate = `total phosphate mg/l`,
    ammonia         = `ammonia mg/l`,
    nitrite         = `nitrite mg/l`,
    nitrate         = `nitrate mg/l`,
    total_nitrogen  = `total nitrogen`,
    chloride        = `chloride mg/l`,
    sulfate         = `sulfate mg/l`
  ) |>
  mutate(system = "Kinzig") |>
  force_numeric()

kinzig_month <- read_march(kinzig_path, "once a month", "date") |>
  rename(
    Cd        = `Cd mg/l`,
    Cu        = `Cu mg/l`,
    Fe        = `Fe mg/l`,
    Ni        = `Ni mg/l`,
    Pb        = `Pb mg/l`,
    Zn        = `Zn mg/l`,
    Ca        = `Ca2+`,
    Mg        = `Mg2+`,
    Na        = `Na+`,
    K         = `K+`,
    carbonate = `Co32-`,
    HCO3      = `HCO3- mg/l`
  ) |>
  mutate(system = "Kinzig") |>
  force_numeric()


# ══════════════════════════════════════════════════════════════════════════════
# 2.  MERGE BOYE + KINZIG
# ══════════════════════════════════════════════════════════════════════════════

all_freq  <- bind_rows(boye_freq,  kinzig_freq)  |> relocate(system, date, site)
all_month <- bind_rows(boye_month, kinzig_month) |> relocate(system, date, site)


# ══════════════════════════════════════════════════════════════════════════════
# 3.  LONG FORMAT + SUMMARY STATS
# ══════════════════════════════════════════════════════════════════════════════

march2025 <- bind_rows(
  to_long(all_freq)  |> mutate(freq = "frequent"),
  to_long(all_month) |> mutate(freq = "monthly")
)

summary_stats <- march2025 |>
  group_by(system, freq, parameter) |>
  summarise(
    n      = n(),
    mean   = round(mean(value),   4),
    median = round(median(value), 4),
    sd     = round(sd(value),     4),
    min    = round(min(value),    4),
    max    = round(max(value),    4),
    .groups = "drop"
  )


# ══════════════════════════════════════════════════════════════════════════════
# 4.  EXCEL OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

wb <- createWorkbook()

sheets <- list(
  "Frequent_Nutrients" = all_freq,
  "Monthly_Metals_Ions" = all_month,
  "Summary_Stats"      = summary_stats,
  "March2025_Long"     = march2025
)

styles <- c("TableStyleMedium9", "TableStyleMedium9",
            "TableStyleMedium2", "TableStyleLight9")

for (i in seq_along(sheets)) {
  sname <- names(sheets)[i]
  addWorksheet(wb, sname)
  writeDataTable(wb, sname, sheets[[i]], tableStyle = styles[i])
  setColWidths(wb, sname, cols = 1:ncol(sheets[[i]]), widths = "auto")
}

saveWorkbook(wb, "March2025_Water_Analysis.xlsx", overwrite = TRUE)
cat("Saved: March2025_Water_Analysis.xlsx\n")





library(readxl)
library(dplyr)
library(openxlsx)

# ══════════════════════════════════════════════════════════════════════════════
# 1.  READ  — confirmed: 1 sheet, 26 cols, duplicate IDs at cols 16-18
# ══════════════════════════════════════════════════════════════════════════════

raw <- read_excel("March2025_Water_Analysis.xlsx", sheet = "Monthly_Metals_Ions")

# Cols 1–15: stream, date, site + metals/ions
metals <- raw[, 1:15]

# Cols 19–26: nutrients only (skip duplicate stream/date/site at 16-18)
nutrients <- raw[, 19:26]
colnames(nutrients) <- c("ortho_phosphate", "total_phosphate", "ammonia",
                         "nitrite", "nitrate", "total_nitrogen",
                         "chloride", "sulfate")

# ══════════════════════════════════════════════════════════════════════════════
# 2.  COMBINE + CLEAN
# ══════════════════════════════════════════════════════════════════════════════

merged <- bind_cols(metals, nutrients) |>
  mutate(date = as.Date(date)) |>
  arrange(stream, date, site) |>
  mutate(across(where(is.numeric), ~ round(., 4)))

# ══════════════════════════════════════════════════════════════════════════════
# 3.  ADD UNITS TO COLUMN NAMES
# ══════════════════════════════════════════════════════════════════════════════

merged <- merged |>
  rename(
    `Cd (mg/l)`              = Cd,
    `Cu (mg/l)`              = Cu,
    `Fe (mg/l)`              = Fe,
    `Ni (mg/l)`              = Ni,
    `Pb (mg/l)`              = Pb,
    `Zn (mg/l)`              = Zn,
    `Ca (mg/l)`              = Ca,
    `Mg (mg/l)`              = Mg,
    `Na (mg/l)`              = Na,
    `K (mg/l)`               = K,
    `carbonate (mg/l)`       = carbonate,
    `HCO3 (mg/l)`            = HCO3,
    `ortho_phosphate (mg/l)` = ortho_phosphate,
    `total_phosphate (mg/l)` = total_phosphate,
    `ammonia (mg/l)`         = ammonia,
    `nitrite (mg/l)`         = nitrite,
    `nitrate (mg/l)`         = nitrate,
    `total_nitrogen (mg/l)`  = total_nitrogen,
    `chloride (mg/l)`        = chloride,
    `sulfate (mg/l)`         = sulfate
  )

# ══════════════════════════════════════════════════════════════════════════════
# 4.  EXCEL OUTPUT — PLAIN, NO TABLE STYLING
# ══════════════════════════════════════════════════════════════════════════════

wb <- createWorkbook()
addWorksheet(wb, "March2025")

writeData(wb, "March2025", merged, startRow = 1, startCol = 1)

addStyle(wb, "March2025",
         createStyle(textDecoration = "bold", border = "Bottom", borderStyle = "medium"),
         rows = 1, cols = 1:ncol(merged), gridExpand = TRUE)

addStyle(wb, "March2025",
         createStyle(numFmt = "YYYY-MM-DD"),
         rows = 2:(nrow(merged) + 1), cols = 2)

setColWidths(wb, "March2025", cols = 1:ncol(merged), widths = "auto")

saveWorkbook(wb, "March2025_Water_Analysis_Clean.xlsx", overwrite = TRUE)
cat("Saved: March2025_Water_Analysis_Clean.xlsx\n")







library(readxl)
library(dplyr)
library(openxlsx)

# ══════════════════════════════════════════════════════════════════════════════
# 1.  READ
# ══════════════════════════════════════════════════════════════════════════════

raw <- read_excel("March2025_Water_Analysis.xlsx", sheet = "Monthly_Metals_Ions")

metals    <- raw[, 1:15]
nutrients <- raw[, 19:26]
colnames(nutrients) <- c("ortho_phosphate", "total_phosphate", "ammonia",
                         "nitrite", "nitrate", "total_nitrogen",
                         "chloride", "sulfate")

merged <- bind_cols(metals, nutrients) |>
  mutate(date = as.Date(date))

# ══════════════════════════════════════════════════════════════════════════════
# 2.  AVERAGE ACROSS THE TWO SAMPLING DATES PER SITE
# ══════════════════════════════════════════════════════════════════════════════

averaged <- merged |>
  group_by(stream, site) |>
  summarise(across(where(is.numeric), ~ round(mean(., na.rm = TRUE), 4)),
            .groups = "drop") |>
  arrange(stream, site)

# ══════════════════════════════════════════════════════════════════════════════
# 3.  ADD UNITS TO COLUMN NAMES
# ══════════════════════════════════════════════════════════════════════════════

averaged <- averaged |>
  rename(
    Cd_mg_l              = Cd,
    Cu_mg_l              = Cu,
    Fe_mg_l              = Fe,
    Ni_mg_l              = Ni,
    Pb_mg_l              = Pb,
    Zn_mg_l              = Zn,
    Ca_mg_l              = Ca,
    Mg_mg_l              = Mg,
    Na_mg_l              = Na,
    K_mg_l               = K,
    carbonate_mg_l       = carbonate,
    HCO3_mg_l            = HCO3,
    ortho_phosphate_mg_l = ortho_phosphate,
    total_phosphate_mg_l = total_phosphate,
    ammonia_mg_l         = ammonia,
    nitrite_mg_l         = nitrite,
    nitrate_mg_l         = nitrate,
    total_nitrogen_mg_l  = total_nitrogen,
    chloride_mg_l        = chloride,
    sulfate_mg_l         = sulfate
  )

# ══════════════════════════════════════════════════════════════════════════════
# 4.  EXCEL OUTPUT — PLAIN
# ══════════════════════════════════════════════════════════════════════════════

wb <- createWorkbook()
addWorksheet(wb, "March2025")

writeData(wb, "March2025", averaged, startRow = 1, startCol = 1)

addStyle(wb, "March2025",
         createStyle(textDecoration = "bold", border = "Bottom", borderStyle = "medium"),
         rows = 1, cols = 1:ncol(averaged), gridExpand = TRUE)

setColWidths(wb, "March2025", cols = 1:ncol(averaged), widths = "auto")

saveWorkbook(wb, "March2025_Water_Analysis_Clean.xlsx", overwrite = TRUE)
cat("Saved: March2025_Water_Analysis_Clean.xlsx\n")