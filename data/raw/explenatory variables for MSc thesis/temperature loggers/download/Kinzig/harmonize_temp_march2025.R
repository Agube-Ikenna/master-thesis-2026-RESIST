# =============================================================================
# Harmonize Water Temperature Data
# Output: one row per stream with Average and SD of Water Temp. [°C]
# NOTE: data in these files runs from May 2021 – October 2024.
#       Change TARGET_MONTH below to any "YYYY-MM" present in the data.
# =============================================================================

library(readxl)
library(dplyr)
library(purrr)
library(lubridate)
library(stringr)

# ── 0. Target month ───────────────────────────────────────────────────────────
TARGET_MONTH <- "2024-10"   # ← change this to any "YYYY-MM" in the data

# ── 1. Set folder path ────────────────────────────────────────────────────────
# Update this path to wherever your Excel files are stored:
folder <- "."   # "." = same directory as this script; or use an absolute path
                # e.g. folder <- "C:/Users/Ikenna/Documents/Kinzig"

# ── 2. Find all .xlsx files (exclude the overview file) ───────────────────────
xlsx_files <- list.files(
  path       = folder,
  pattern    = "\\.xlsx$",
  full.names = TRUE
)
xlsx_files <- xlsx_files[!grepl("_OVERVIEW", xlsx_files)]

cat("Found", length(xlsx_files), "data files.\n")

# ── 3. Helper: read one sheet, return tidy tibble ─────────────────────────────
read_logger_sheet <- function(file, sheet) {

  df <- tryCatch(
    read_excel(file, sheet = sheet, col_types = "text"),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0)      return(NULL)
  # Find the temperature column (header has a literal newline: "Water Temp.\n[°C]")
  temp_col <- grep("Water Temp", names(df), value = TRUE, fixed = FALSE)[1]
  if (is.na(temp_col))          return(NULL)
  if (!"date" %in% names(df))   return(NULL)

  stream   <- str_remove(sheet, "_[AB]$")
  raw_date <- df[["date"]]
  raw_temp <- df[[temp_col]]

  # Parse date → plain character "YYYY-MM-DD" (never a Date/POSIXct object).
  # Keeping it as character means bind_rows() in map_dfr() has nothing to
  # coerce and cannot create a list-column.
  #   Excel serial stored as text (e.g. "45717") → numeric origin conversion
  #   ISO / datetime string (e.g. "2025-03-01 00:00:00") → first 10 chars
  num_serial  <- suppressWarnings(as.numeric(raw_date))
  date_char   <- ifelse(
    !is.na(num_serial),
    as.character(as.Date(num_serial, origin = "1899-12-30")),
    substr(raw_date, 1, 10)
  )  # result: "YYYY-MM-DD" or NA — always a plain character vector

  out <- data.frame(
    stream = stream,
    date   = date_char,
    temp   = suppressWarnings(as.numeric(raw_temp)),
    stringsAsFactors = FALSE
  )

  out[!is.na(out$date) & nchar(out$date) == 10 & !is.na(out$temp), ]
}

# ── 4. Loop over all files and all data sheets ────────────────────────────────
all_data <- map_dfr(xlsx_files, function(file) {

  sheets <- excel_sheets(file)
  # Keep only logger sheets (skip "protocol")
  data_sheets <- sheets[sheets != "protocol"]

  map_dfr(data_sheets, ~ read_logger_sheet(file, .x))
})

cat("Total rows loaded:", nrow(all_data), "\n")

# ── 5. Filter to target month ─────────────────────────────────────────────────
# date is a plain character "YYYY-MM-DD" — startsWith() needs no type conversion
month_data <- all_data[startsWith(all_data$date, TARGET_MONTH), ]

cat("Rows in", TARGET_MONTH, ":", nrow(month_data), "\n")

# ── 6. Summarise: Average and SD per stream ───────────────────────────────────
result <- month_data %>%
  group_by(Stream = stream) %>%
  summarise(
    Date    = TARGET_MONTH,
    Average = round(mean(temp, na.rm = TRUE), 3),
    SD      = round(sd(temp,   na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(Stream)

# ── 7. Preview ────────────────────────────────────────────────────────────────
print(result)

# ── 8. Export (optional — uncomment to save) ──────────────────────────────────
# write.csv(result, file.path(folder, "WaterTemp_March2025_Summary.csv"),
#           row.names = FALSE)
