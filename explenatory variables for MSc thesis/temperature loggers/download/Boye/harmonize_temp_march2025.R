# =============================================================================
# Harmonize Water Temperature Logger Data – Monthly Summary (March 2025)
# =============================================================================
# Required packages: readxl, dplyr, lubridate, purrr, stringr
# Install if needed:
#   install.packages(c("readxl", "dplyr", "lubridate", "purrr", "stringr"))
# =============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(purrr)
library(stringr)

# ── 1. Locate the folder ──────────────────────────────────────────────────────
# Option A: set path explicitly (edit the string below to match your machine)
folder <- ""   # e.g. "C:/Users/ikenna/Documents/Boye"

# Option B: if folder is left blank above, a dialog box will open so you can
#           browse to the folder interactively (requires RStudio or tcltk)
if (!nzchar(folder) || !dir.exists(folder)) {
  message("Folder not set or not found – opening folder picker…")
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    folder <- rstudioapi::selectDirectory(caption = "Select the Boye data folder")
  } else {
    folder <- tcltk::tk_choose.dir(caption = "Select the Boye data folder")
  }
}

if (is.null(folder) || !nzchar(folder) || !dir.exists(folder)) {
  stop("No valid folder selected. Please set the 'folder' variable at the top of the script.")
}

cat("Using folder:", folder, "\n")

# ── 2. Collect all Excel files (excluding the OVERVIEW file) ─────────────────
files <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)
files <- files[!grepl("OVERVIEW", files, ignore.case = TRUE)]

if (length(files) == 0) {
  stop("No .xlsx files found in: ", folder,
       "\nCheck that the folder path is correct.")
}
cat("Found", length(files), "data files.\n")

# ── 3. Helper: read one sheet by column position ──────────────────────────────
# col_names=FALSE + skip=1 skips the header row entirely.
# Column layout: 1=LoggerName, 2=LoggerID, 3=date, 4=year, 5=month,
#                6=time, 7=Water Temp [°C]
read_logger_sheet <- function(file, sheet) {

  raw <- read_excel(file, sheet = sheet,
                    col_names = FALSE,
                    skip = 1)

  if (ncol(raw) < 7) {
    stop("Sheet has fewer than 7 columns (found ", ncol(raw), ")")
  }

  tibble(
    stream       = str_remove(sheet, "_[AB]$"),
    date         = as.Date(raw[[3]]),
    water_temp_C = suppressWarnings(as.numeric(raw[[7]]))
  )
}

# ── 4. Read every data sheet from every file ─────────────────────────────────
all_data <- map_dfr(files, function(file) {

  sheets <- excel_sheets(file)
  sheets <- sheets[sheets != "protocol"]

  map_dfr(sheets, function(sheet) {
    tryCatch(
      read_logger_sheet(file, sheet),
      error = function(e) {
        message("  Skipped ", basename(file), " / ", sheet, ": ", e$message)
        NULL
      }
    )
  })
})

cat("Total rows read:", nrow(all_data), "\n")
cat("Columns in all_data:", paste(names(all_data), collapse = ", "), "\n")

if (nrow(all_data) == 0 || !"water_temp_C" %in% names(all_data)) {
  stop("No data was loaded. Check the 'Skipped' messages above for errors per sheet.")
}

# ── 5. Filter for March 2025 & compute monthly statistics ────────────────────
result <- all_data %>%
  filter(!is.na(water_temp_C), !is.na(date)) %>%
  mutate(
    yr      = year(date),
    mo      = month(date),
    date_ym = format(date, "%Y-%m")
  ) %>%
  filter(yr == 2025, mo == 3) %>%
  group_by(Stream = stream, Date = date_ym) %>%
  summarise(
    Average = round(mean(water_temp_C, na.rm = TRUE), 3),
    SD      = round(sd(water_temp_C,   na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(Stream)

cat("\nRows matching March 2025:", nrow(result), "\n\n")

# ── 6. Print result ───────────────────────────────────────────────────────────
print(result, n = Inf)

# ── 7. (Optional) Export to CSV ───────────────────────────────────────────────
# Uncomment to save:
 write.csv(result, file.path(folder, "monthly_summary_march2025.csv"), row.names = FALSE)
