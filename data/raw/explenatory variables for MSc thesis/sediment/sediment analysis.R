library(readxl)
library(dplyr)
library(writexl)

# ── Load ──────────────────────────────────────────────────────────────────────
raw <- read_excel(
  "SFB_FieldStudies_SubstrateEstimation_Boye_2026-01-22 (inkl. sampling data 2025).xlsx",
  sheet = "year 2025",
  col_names = FALSE
)

# ── Remove fully empty rows and columns ───────────────────────────────────────
raw <- raw[rowSums(!is.na(raw)) > 0, ]
raw <- raw[, colSums(!is.na(raw)) > 0]

# ── Assign column names from two header rows ──────────────────────────────────
col_names_r1 <- as.character(raw[1, ])
col_names_r2 <- as.character(raw[2, ])
col_names <- ifelse(!is.na(col_names_r2) & col_names_r2 != "NA", col_names_r2, col_names_r1)
col_names[is.na(col_names) | col_names == "NA"] <- paste0("V", which(is.na(col_names) | col_names == "NA"))
colnames(raw) <- col_names

# ── Keep data rows only ───────────────────────────────────────────────────────
df <- raw[3:nrow(raw), ]
df <- df[!((!is.na(df$`yyyy-mm-dd`)) & grepl("x = share", df$`yyyy-mm-dd`, ignore.case = TRUE)), ]

# ── Mark unsampled site (before mega is removed) ─────────────────────────────
not_sampled <- !is.na(df$mega) & grepl("not sampled", df$mega, ignore.case = TRUE)

# ── Convert "x" to 0, then to numeric (on all substrate cols at this stage) ──
substrate_cols_raw <- setdiff(
  colnames(df),
  c("sites", "yyyy-mm-dd", "sum", "notes on habitats", "notes on biology")
)

df[substrate_cols_raw] <- lapply(df[substrate_cols_raw], function(col) {
  col <- trimws(as.character(col))
  col[col == "x"] <- "0"
  suppressWarnings(as.numeric(col))
})

# Set unsampled site to NA
df[not_sampled, substrate_cols_raw] <- NA

# ── Convert date and sum ──────────────────────────────────────────────────────
df$`yyyy-mm-dd` <- as.Date(as.numeric(df$`yyyy-mm-dd`), origin = "1899-12-30")
df$sum          <- as.numeric(df$sum)

# ── Drop all unwanted columns in one step ─────────────────────────────────────
drop_cols <- c("mega", "arg", "tech_2", "algae", "sapro", "debris",
               "notes on habitats", "notes on biology")
df <- df %>% select(-any_of(drop_cols))

# ── Define final substrate cols AFTER dropping ────────────────────────────────
substrate_cols <- setdiff(colnames(df), c("sites", "yyyy-mm-dd", "sum"))

# Replace remaining NA with 0 for sampled sites only
sampled_rows <- !not_sampled

df[sampled_rows, substrate_cols] <- lapply(
  df[sampled_rows, substrate_cols],
  function(col) replace(col, is.na(col), 0)
)

# ── Save clean data ───────────────────────────────────────────────────────────
write_xlsx(df, "Boye_2025_clean.xlsx")

# ── Summary statistics ────────────────────────────────────────────────────────
summary_stats <- data.frame(
  variable  = substrate_cols,
  n         = sapply(df[substrate_cols], function(x) sum(!is.na(x))),
  n_zero    = sapply(df[substrate_cols], function(x) sum(x == 0,  na.rm = TRUE)),
  n_present = sapply(df[substrate_cols], function(x) sum(x > 0,   na.rm = TRUE)),
  mean      = sapply(df[substrate_cols], function(x) round(mean(x,   na.rm = TRUE), 2)),
  median    = sapply(df[substrate_cols], function(x) round(median(x, na.rm = TRUE), 2)),
  sd        = sapply(df[substrate_cols], function(x) round(sd(x,     na.rm = TRUE), 2)),
  min       = sapply(df[substrate_cols], function(x) round(min(x,    na.rm = TRUE), 2)),
  max       = sapply(df[substrate_cols], function(x) round(max(x,    na.rm = TRUE), 2)),
  row.names = NULL
)

write_xlsx(summary_stats, "Boye_2025_summary_stats.xlsx")

# ── Print ─────────────────────────────────────────────────────────────────────
print(as.data.frame(df), row.names = FALSE)
print(summary_stats)







library(readxl)
library(dplyr)
library(writexl)

# ── Load ──────────────────────────────────────────────────────────────────────
raw <- read_excel(
  "SFB_FieldStudies_SubstrateEstimation_Kinzig_2025-05-12 (inkl. sampling date 2025).xlsx",
  sheet = "year 2025",
  col_names = FALSE
)

# ── Remove fully empty rows and columns ───────────────────────────────────────
raw <- raw[rowSums(!is.na(raw)) > 0, ]
raw <- raw[, colSums(!is.na(raw)) > 0]

# ── Assign column names from two header rows ──────────────────────────────────
col_names_r1 <- as.character(raw[1, ])
col_names_r2 <- as.character(raw[2, ])
col_names <- ifelse(!is.na(col_names_r2) & col_names_r2 != "NA", col_names_r2, col_names_r1)
col_names[is.na(col_names) | col_names == "NA"] <- paste0("V", which(is.na(col_names) | col_names == "NA"))
colnames(raw) <- col_names

# ── Keep data rows only ───────────────────────────────────────────────────────
df <- raw[3:nrow(raw), ]

# Drop legend row ("x = share below 5%")
df <- df[!((!is.na(df$`yyyy-mm-dd`)) & grepl("x = share", df$`yyyy-mm-dd`, ignore.case = TRUE)), ]

# Drop trailing note rows (rows without a site code)
df <- df[!is.na(df$sites), ]

# ── Mark unsampled sites (none expected in Kinzig but kept for robustness) ────
not_sampled <- !is.na(df$mega) & grepl("not sampled", df$mega, ignore.case = TRUE)

# ── Convert "-" and "x" to 0, then to numeric ────────────────────────────────
substrate_cols_raw <- setdiff(
  colnames(df),
  c("sites", "yyyy-mm-dd", "sum", "notes on habitats", "notes on biology")
)

df[substrate_cols_raw] <- lapply(df[substrate_cols_raw], function(col) {
  col <- trimws(as.character(col))
  col[col == "x"] <- "0"   # trace presence (<5%)
  col[col == "-"] <- "0"   # absent substrate
  suppressWarnings(as.numeric(col))
})

# Set unsampled site rows to NA (precautionary)
df[not_sampled, substrate_cols_raw] <- NA

# ── Convert date and sum ──────────────────────────────────────────────────────
df$`yyyy-mm-dd` <- as.Date(as.numeric(df$`yyyy-mm-dd`), origin = "1899-12-30")
df$sum          <- as.numeric(df$sum)

# ── Drop all unwanted columns in one step ─────────────────────────────────────
drop_cols <- c("mega", "arg", "tech_2", "algae", "sapro", "debris",
               "notes on habitats", "notes on biology")
df <- df %>% select(-any_of(drop_cols))

# ── Define final substrate cols AFTER dropping ────────────────────────────────
substrate_cols <- setdiff(colnames(df), c("sites", "yyyy-mm-dd", "sum"))

# ── Replace remaining NA with 0 for sampled sites only ───────────────────────
sampled_rows <- !not_sampled
df[sampled_rows, substrate_cols] <- lapply(
  df[sampled_rows, substrate_cols],
  function(col) replace(col, is.na(col), 0)
)

# ── Save clean data ───────────────────────────────────────────────────────────
write_xlsx(df, "Kinzig_2025_clean.xlsx")

# ── Summary statistics ────────────────────────────────────────────────────────
summary_stats <- data.frame(
  variable  = substrate_cols,
  n         = sapply(df[substrate_cols], function(x) sum(!is.na(x))),
  n_zero    = sapply(df[substrate_cols], function(x) sum(x == 0,  na.rm = TRUE)),
  n_present = sapply(df[substrate_cols], function(x) sum(x > 0,   na.rm = TRUE)),
  mean      = sapply(df[substrate_cols], function(x) round(mean(x,   na.rm = TRUE), 2)),
  median    = sapply(df[substrate_cols], function(x) round(median(x, na.rm = TRUE), 2)),
  sd        = sapply(df[substrate_cols], function(x) round(sd(x,     na.rm = TRUE), 2)),
  min       = sapply(df[substrate_cols], function(x) round(min(x,    na.rm = TRUE), 2)),
  max       = sapply(df[substrate_cols], function(x) round(max(x,    na.rm = TRUE), 2)),
  row.names = NULL
)

write_xlsx(summary_stats, "Kinzig_2025_summary_stats.xlsx")

# ── Print ─────────────────────────────────────────────────────────────────────
print(as.data.frame(df), row.names = FALSE)
print(summary_stats)



















