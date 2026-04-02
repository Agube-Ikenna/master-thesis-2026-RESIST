library(readxl)
library(writexl)
library(dplyr)
library(tidyr)

# Load raw data (skip title row, use row 2 as header)
df_raw <- read_excel(
  "SFB_MACROINVERTEBRATE_BOYE_2O25.xlsx",
  sheet = "Sheet1",
  col_names = FALSE
)

# Extract stream names and dates (row 2 = headers, row 3 = dates)
stream_names <- as.character(df_raw[2, 5:ncol(df_raw)])
stream_names <- gsub("\n", " ", stream_names)  # clean line breaks

dates <- as.Date(as.numeric(df_raw[3, 5:ncol(df_raw)]), origin = "1899-12-30")

# Extract taxon data (rows 4 onwards)
data_block <- df_raw[4:nrow(df_raw), ]

taxon_names <- as.character(data_block[[4]])  # Column D = Taxon name
values      <- data_block[, 5:ncol(data_block)]  # Stream value columns

# Convert values to numeric matrix
values_mat <- apply(values, 2, as.numeric)
values_mat[is.na(values_mat)] <- 0

# Transpose: streams on rows, taxa on columns
transposed <- as.data.frame(t(values_mat))
colnames(transposed) <- taxon_names

# Add stream name and date as first two columns
result <- bind_cols(
  data.frame(Stream = stream_names, Date_2025 = dates),
  transposed
)

# Write output
write_xlsx(result, "Streams_on_Rows_2025.xlsx")




library(readxl)
library(writexl)
library(dplyr)

# Load raw data
df_raw <- read_excel(
  "SFB_MACROINVERTEBRATE_BOYE_2O25.xlsx",
  sheet = "Sheet1",
  col_names = FALSE
)

# Extract components
stream_names <- gsub("\n", " ", as.character(df_raw[2, 5:ncol(df_raw)]))
dates        <- as.Date(as.numeric(df_raw[3, 5:ncol(df_raw)]), origin = "1899-12-30")
taxon_names  <- as.character(df_raw[4:nrow(df_raw), 4, drop = TRUE])
values_mat   <- apply(df_raw[4:nrow(df_raw), 5:ncol(df_raw)], 2, as.numeric)
values_mat[is.na(values_mat)] <- 0

# Transpose: streams on rows, taxa on columns
transposed <- as.data.frame(t(values_mat))
colnames(transposed) <- taxon_names

# Combine stream name + date + values
result <- bind_cols(
  data.frame(Stream = stream_names, Date_2025 = dates, check.names = FALSE),
  transposed
)

# Save
write_xlsx(result, "Streams_on_Rows_2025.xlsx")

cat("Done:", nrow(result), "streams x", ncol(result), "columns\n")








library(readxl)
library(writexl)
library(dplyr)

# Load raw data — treat "-" and blanks as NA, then replace with 0
df_raw <- read_excel(
  "SFB_MACROINVERTEBRATE_BOYE_2O25.xlsx",
  sheet    = "Sheet1",
  col_names = FALSE,
  na       = c("", "NA", "-")
)

# Extract components
stream_names <- gsub("\n", " ", as.character(df_raw[2, 5:ncol(df_raw)]))
dates        <- as.Date(as.numeric(df_raw[3, 5:ncol(df_raw)]), origin = "1899-12-30")
taxon_names  <- as.character(df_raw[4:nrow(df_raw), 4, drop = TRUE])
values_mat   <- apply(df_raw[4:nrow(df_raw), 5:ncol(df_raw)], 2, as.numeric)

# Replace all NA / "-" with 0
values_mat[is.na(values_mat)] <- 0

# Transpose: streams on rows, taxa on columns
transposed <- as.data.frame(t(values_mat))
colnames(transposed) <- taxon_names

# Combine stream name + date + values
result <- bind_cols(
  data.frame(Stream = stream_names, Date_2025 = dates, check.names = FALSE),
  transposed
)

# Save
write_xlsx(result, "Streams_on_Rows_2025.xlsx")

cat("Done:", nrow(result), "streams x", ncol(result) - 2, "taxa columns\n")







library(readxl)
library(writexl)
library(dplyr)

# Load raw data
df_raw <- read_excel(
  "SFB_MACROINVERTEBRATE_BOYE_2O25.xlsx",
  sheet     = "Sheet1",
  col_names = FALSE,
  na        = c("", "NA", "-")
)

# Extract components
stream_names <- gsub("\n", " ", as.character(df_raw[2, 5:ncol(df_raw)]))
dates        <- as.Date(as.numeric(df_raw[3, 5:ncol(df_raw)]), origin = "1899-12-30")
taxon_names  <- as.character(df_raw[4:nrow(df_raw), 4, drop = TRUE])

# Read values as-is, replace NA with 0, no numeric coercion that adds digits
values_block <- df_raw[4:nrow(df_raw), 5:ncol(df_raw)]
values_block[is.na(values_block)] <- 0

# Convert preserving original precision
values_mat <- as.matrix(sapply(values_block, function(x) {
  v <- suppressWarnings(as.numeric(x))
  v[is.na(v)] <- 0
  v
}))

# Transpose: streams on rows, taxa on columns
transposed <- as.data.frame(t(values_mat), check.names = FALSE)
colnames(transposed) <- taxon_names

# Combine
result <- bind_cols(
  data.frame(Stream = stream_names, Date_2025 = dates, check.names = FALSE),
  transposed
)

# Save — numbers written exactly as read, no formatting added
write_xlsx(result, "Streams_on_Rows_2025.xlsx")

cat("Done:", nrow(result), "streams x", ncol(result) - 2, "taxa columns\n")


library(readxl)
library(writexl)
library(dplyr)

# ── Load raw data ──────────────────────────────────────────────────────────────
df_raw <- read_excel(
  "SFB_MACROINVERTEBRATE_BOYE_2O25.xlsx",
  sheet     = "Sheet1",
  col_names = FALSE,
  na        = c("", "NA", "-")
)

stream_names <- gsub("\n", " ", as.character(df_raw[2, 5:ncol(df_raw)]))
dates        <- as.Date(as.numeric(df_raw[3, 5:ncol(df_raw)]), origin = "1899-12-30")
family       <- as.character(df_raw[4:nrow(df_raw), 2, drop = TRUE])
taxon_name   <- as.character(df_raw[4:nrow(df_raw), 4, drop = TRUE])

values_block <- df_raw[4:nrow(df_raw), 5:ncol(df_raw)]
values_mat   <- matrix(
  as.numeric(unlist(values_block)),
  nrow = nrow(values_block),
  ncol = ncol(values_block)
)
values_mat[is.na(values_mat)] <- 0
# values_mat: rows = taxa, cols = streams

n_streams <- length(stream_names)
n_taxa    <- length(taxon_name)

# ── Per-stream statistics ──────────────────────────────────────────────────────
stream_totals  <- colSums(values_mat)
stream_mean    <- colMeans(values_mat)
stream_median  <- apply(values_mat, 2, median)

# Highest species per stream
highest_sp_idx <- apply(values_mat, 2, which.max)
highest_species <- taxon_name[highest_sp_idx]
highest_sp_val  <- values_mat[cbind(highest_sp_idx, seq_len(n_streams))]

# Highest family per stream
top_family_for_stream <- function(stream_idx) {
  fam_totals <- tapply(values_mat[, stream_idx], family, sum, na.rm = TRUE)
  fam_totals <- fam_totals[names(fam_totals) != "NA"]
  top        <- names(which.max(fam_totals))
  list(name = top, value = fam_totals[[top]])
}

top_families    <- lapply(seq_len(n_streams), top_family_for_stream)
highest_family  <- sapply(top_families, `[[`, "name")
highest_fam_val <- sapply(top_families, `[[`, "value")

# ── Sheet 1: Per-Stream Summary ────────────────────────────────────────────────
sheet1 <- data.frame(
  Stream                        = stream_names,
  Sampling_Date                 = as.character(dates),
  Total_Abundance               = stream_totals,
  Mean_Abundance_per_taxon      = stream_mean,
  Median_Abundance_per_taxon    = stream_median,
  Highest_Family                = highest_family,
  Family_Total_Abundance        = highest_fam_val,
  Highest_Species               = highest_species,
  Species_Abundance             = highest_sp_val,
  check.names = FALSE,
  row.names   = NULL
)

# ── Sheet 2: Catchment Total ───────────────────────────────────────────────────
catchment_total  <- sum(values_mat)
catchment_mean   <- mean(values_mat)
catchment_median <- median(values_mat)

taxa_totals_all  <- rowSums(values_mat)
top_sp           <- taxon_name[which.max(taxa_totals_all)]
top_sp_val       <- max(taxa_totals_all)

fam_totals_all <- tapply(taxa_totals_all, family, sum, na.rm = TRUE)
fam_totals_all <- fam_totals_all[names(fam_totals_all) != "NA"]
top_fam        <- names(which.max(fam_totals_all))
top_fam_val    <- max(fam_totals_all)

sheet2 <- data.frame(
  Metric = c(
    "Total Macroinvertebrate Abundance (all streams)",
    "Mean Abundance per taxon per stream",
    "Median Abundance per taxon per stream",
    "Dominant Family (highest total)",
    "Dominant Species (highest total)",
    "Number of Streams sampled",
    "Number of Taxa recorded"
  ),
  Value = c(
    catchment_total, catchment_mean, catchment_median,
    top_fam, top_sp,
    n_streams, n_taxa
  ),
  Notes = c(
    paste("Sum across all", n_streams, "streams &", n_taxa, "taxa"),
    paste("Average of all", n_taxa, "x", n_streams, "=", n_taxa * n_streams, "values"),
    paste("Median of all",  n_taxa, "x", n_streams, "=", n_taxa * n_streams, "values"),
    paste("Total:", round(top_fam_val, 3)),
    paste("Total:", round(top_sp_val,  3)),
    "Year 2025",
    "Across all streams"
  ),
  check.names = FALSE
)

# ── Sheet 3: Top 10 Species & Families ────────────────────────────────────────
top10_species <- data.frame(
  Rank            = 1:10,
  Species         = taxon_name[order(taxa_totals_all, decreasing = TRUE)[1:10]],
  Total_Abundance = sort(taxa_totals_all, decreasing = TRUE)[1:10],
  check.names = FALSE, row.names = NULL
)

top10_families <- data.frame(
  Rank            = 1:10,
  Family          = names(sort(fam_totals_all, decreasing = TRUE))[1:10],
  Total_Abundance = sort(fam_totals_all, decreasing = TRUE)[1:10],
  check.names = FALSE, row.names = NULL
)

# ── Save all sheets ────────────────────────────────────────────────────────────
write_xlsx(
  list(
    "Per-Stream Summary"    = sheet1,
    "Catchment Total"       = sheet2,
    "Top Species"           = top10_species,
    "Top Families"          = top10_families
  ),
  "Macroinvertebrate_Summary_2025.xlsx"
)

cat("Done!\n")
cat("Sheet 1:", nrow(sheet1), "streams\n")
cat("Catchment total abundance:", round(catchment_total, 3), "\n")
cat("Dominant species:", top_sp, "\n")
cat("Dominant family:", top_fam, "\n")


