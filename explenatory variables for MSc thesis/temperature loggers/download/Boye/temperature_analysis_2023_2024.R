# =============================================================================
#  Boye Catchment – Water Temperature Analysis  |  2023 & 2024
# =============================================================================
#  What this script produces
#  ─────────────────────────
#  STATISTICS (data frames)
#    daily_stats    – daily mean / min / max / range per stream
#    monthly_stats  – monthly mean / min / max / range / SD per stream
#    yearly_stats   – yearly  mean / min / max / range / SD per stream
#    catch_daily    – same metrics aggregated across ALL streams (catchment)
#    catch_monthly  – catchment monthly summary
#    catch_yearly   – catchment yearly summary
#    metadata       – data-coverage and sensor info per stream
#
#  CHARTS (saved to sub-folder  <folder>/plots/)
#    01_monthly_heatmap        – stream × month mean-temp heatmap
#    02_monthly_lines          – monthly averages per stream (Jan–Dec)
#    03_catchment_trend        – catchment monthly avg + linear & exponential fits
#    04_daily_range_boxplots   – daily temperature range distribution by month
#    05_yearly_comparison      – 2023 vs 2024 mean temp per stream
#    06_stream_trends          – individual stream monthly trends with fits
#    07_catchment_annual_cycle – catchment daily mean ± SD ribbon
#
#  OPTIONAL EXCEL EXPORT (Section 9)
#    One .xlsx with sheets: daily_stats, monthly_stats, yearly_stats,
#    catchment_daily, catchment_monthly, catchment_yearly, metadata
#
#  Required packages:
#    install.packages(c("readxl","dplyr","lubridate","purrr","stringr",
#                       "ggplot2","scales","viridis","openxlsx","tidyr","broom"))
# =============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(purrr)
library(stringr)
library(ggplot2)
library(scales)
library(viridis)
library(tidyr)
library(broom)

# ── 0. Folder paths ───────────────────────────────────────────────────────────
# Set to the folder containing your Excel files
folder <- ""   # e.g. "C:/Users/ikenna/Documents/Boye"

# ── Auto-detect if left blank ──────────────────────────────────────────────────
if (!nzchar(folder) || !dir.exists(folder)) {
  message("No folder set – opening folder picker…")
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    folder <- rstudioapi::selectDirectory(caption = "Select the Boye data folder")
  } else {
    folder <- tcltk::tk_choose.dir(caption = "Select the Boye data folder")
  }
}
stopifnot("Folder not found – set the 'folder' variable." = dir.exists(folder))

plots_dir <- file.path(folder, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# =============================================================================
# SECTION 1 – Load & clean raw data (2023–2024 only)
# =============================================================================
files <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)
files <- files[!grepl("OVERVIEW|analysis", basename(files), ignore.case = TRUE)]
cat("Loading", length(files), "files…\n")

read_sheet <- function(file, sheet) {
  raw <- read_excel(file, sheet = sheet, col_names = FALSE, skip = 1)
  if (ncol(raw) < 7) return(NULL)
  tibble(
    stream       = str_remove(sheet, "_[AB]$"),
    date         = as.Date(raw[[3]]),
    water_temp_C = suppressWarnings(as.numeric(raw[[7]]))
  )
}

raw_all <- map_dfr(files, function(f) {
  sheets <- excel_sheets(f)
  sheets <- sheets[sheets != "protocol"]
  map_dfr(sheets, ~ tryCatch(read_sheet(f, .x),
                             error = function(e) { message("Skip: ", .x); NULL }))
}) %>%
  filter(!is.na(water_temp_C), !is.na(date)) %>%
  mutate(
    year  = year(date),
    month = month(date),
    month_label = factor(month.abb[month], levels = month.abb),  # Jan–Dec order
    day   = yday(date)   # day of year
  ) %>%
  filter(year %in% c(2023, 2024))   # ← 2023 and 2024 only

cat("Rows loaded (2023–2024):", nrow(raw_all), "\n")
cat("Streams found:", n_distinct(raw_all$stream), "\n")

# =============================================================================
# SECTION 2 – Daily statistics per stream
# =============================================================================
daily_stats <- raw_all %>%
  group_by(stream, year, date, month, month_label) %>%
  summarise(
    daily_mean  = round(mean(water_temp_C,  na.rm = TRUE), 3),
    daily_min   = round(min(water_temp_C,   na.rm = TRUE), 3),
    daily_max   = round(max(water_temp_C,   na.rm = TRUE), 3),
    daily_range = round(daily_max - daily_min, 3),
    n_readings  = n(),
    .groups = "drop"
  ) %>%
  arrange(stream, date)

cat("Daily stats rows:", nrow(daily_stats), "\n")

# =============================================================================
# SECTION 3 – Monthly statistics per stream
# =============================================================================
monthly_stats <- daily_stats %>%
  group_by(stream, year, month, month_label) %>%
  summarise(
    monthly_mean       = round(mean(daily_mean,  na.rm = TRUE), 3),
    monthly_min        = round(min(daily_min,    na.rm = TRUE), 3),
    monthly_max        = round(max(daily_max,    na.rm = TRUE), 3),
    monthly_range      = round(monthly_max - monthly_min, 3),
    monthly_sd         = round(sd(daily_mean,    na.rm = TRUE), 3),
    monthly_mean_range = round(mean(daily_range, na.rm = TRUE), 3),  # avg daily range
    n_days             = n(),
    .groups = "drop"
  ) %>%
  arrange(stream, year, month)

# =============================================================================
# SECTION 4 – Yearly statistics per stream
# =============================================================================
yearly_stats <- monthly_stats %>%
  group_by(stream, year) %>%
  summarise(
    yearly_mean       = round(mean(monthly_mean,  na.rm = TRUE), 3),
    yearly_min        = round(min(monthly_min,    na.rm = TRUE), 3),
    yearly_max        = round(max(monthly_max,    na.rm = TRUE), 3),
    yearly_range      = round(yearly_max - yearly_min, 3),
    yearly_sd         = round(sd(monthly_mean,    na.rm = TRUE), 3),
    yearly_mean_range = round(mean(monthly_range, na.rm = TRUE), 3),
    months_covered    = n(),
    .groups = "drop"
  ) %>%
  arrange(stream, year)

# =============================================================================
# SECTION 5 – Catchment-level aggregations  (all streams combined)
# =============================================================================
catch_daily <- raw_all %>%
  group_by(year, date, month, month_label) %>%
  summarise(
    catch_mean  = round(mean(water_temp_C, na.rm = TRUE), 3),
    catch_min   = round(min(water_temp_C,  na.rm = TRUE), 3),
    catch_max   = round(max(water_temp_C,  na.rm = TRUE), 3),
    catch_range = round(catch_max - catch_min, 3),
    n_streams   = n_distinct(stream),
    .groups = "drop"
  ) %>%
  arrange(date)

catch_monthly <- catch_daily %>%
  group_by(year, month, month_label) %>%
  summarise(
    catch_monthly_mean  = round(mean(catch_mean,  na.rm = TRUE), 3),
    catch_monthly_min   = round(min(catch_min,    na.rm = TRUE), 3),
    catch_monthly_max   = round(max(catch_max,    na.rm = TRUE), 3),
    catch_monthly_range = round(catch_monthly_max - catch_monthly_min, 3),
    catch_monthly_sd    = round(sd(catch_mean,    na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(year, month)

catch_yearly <- catch_monthly %>%
  group_by(year) %>%
  summarise(
    catch_yearly_mean  = round(mean(catch_monthly_mean, na.rm = TRUE), 3),
    catch_yearly_min   = round(min(catch_monthly_min,   na.rm = TRUE), 3),
    catch_yearly_max   = round(max(catch_monthly_max,   na.rm = TRUE), 3),
    catch_yearly_range = round(catch_yearly_max - catch_yearly_min, 3),
    catch_yearly_sd    = round(sd(catch_monthly_mean,   na.rm = TRUE), 3),
    .groups = "drop"
  )

cat("\n── Catchment Yearly Summary ──\n")
print(catch_yearly)

# =============================================================================
# SECTION 6 – Metadata
# =============================================================================
metadata <- raw_all %>%
  group_by(stream) %>%
  summarise(
    years_covered    = paste(sort(unique(year)), collapse = ", "),
    date_start       = min(date),
    date_end         = max(date),
    total_readings   = n(),
    total_days       = n_distinct(date),
    mean_temp_all    = round(mean(water_temp_C, na.rm = TRUE), 3),
    min_temp_all     = round(min(water_temp_C,  na.rm = TRUE), 3),
    max_temp_all     = round(max(water_temp_C,  na.rm = TRUE), 3),
    overall_range    = round(max_temp_all - min_temp_all, 3),
    pct_missing_days = round(100 * (1 - total_days / as.numeric(date_end - date_start + 1)), 1),
    .groups = "drop"
  ) %>%
  arrange(stream)

cat("\n── Metadata preview ──\n")
print(metadata, n = Inf)

# =============================================================================
# SECTION 7 – Visualisations
# =============================================================================
# Shared theme
theme_boye <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    legend.position  = "right",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

month_pal <- viridis(12, option = "C")  # 12 colours for months

# ── 7.1 Monthly heatmap: stream × month ───────────────────────────────────────
p1 <- monthly_stats %>%
  ggplot(aes(x = month_label, y = stream, fill = monthly_mean)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", monthly_mean)), size = 2.5, colour = "white") +
  facet_wrap(~ year, ncol = 1) +
  scale_fill_viridis_c(option = "C", name = "Mean\nTemp (°C)") +
  labs(title   = "Monthly Mean Water Temperature – Boye Catchment",
       subtitle = "Each cell = monthly average (°C) | faceted by year",
       x = NULL, y = NULL) +
  theme_boye +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(plots_dir, "01_monthly_heatmap.png"), p1,
       width = 14, height = 10, dpi = 180)
cat("Saved plot 1\n")

# ── 7.2 Monthly line chart: all streams Jan–Dec ───────────────────────────────
p2 <- monthly_stats %>%
  ggplot(aes(x = month_label, y = monthly_mean,
             group = interaction(stream, year), colour = stream)) +
  geom_line(alpha = 0.7, linewidth = 0.7) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ year, ncol = 1) +
  scale_colour_viridis_d(option = "D", name = "Stream") +
  labs(title   = "Monthly Average Temperature per Stream (Jan–Dec)",
       subtitle = "2023 and 2024 | Boye Catchment",
       x = NULL, y = "Mean Temperature (°C)") +
  theme_boye

ggsave(file.path(plots_dir, "02_monthly_lines.png"), p2,
       width = 14, height = 9, dpi = 180)
cat("Saved plot 2\n")

# ── 7.3 Catchment trend with linear & exponential fits ────────────────────────
# Prepare data: one row per year–month for catchment
catch_plot <- catch_monthly %>%
  mutate(month_num = month,
         year_f    = factor(year))

# Fit models per year
fit_linear <- function(df) {
  lm(catch_monthly_mean ~ month_num, data = df)
}
fit_exponen <- function(df) {
  # log-linear: log(T) = a + b*month  →  T = e^(a+b*m)
  # Shift temp to ensure all values > 0 before log
  df <- df %>% mutate(temp_pos = catch_monthly_mean - min(catch_monthly_mean) + 0.1)
  tryCatch(lm(log(temp_pos) ~ month_num, data = df), error = function(e) NULL)
}

pred_data <- catch_plot %>%
  group_by(year) %>%
  group_map(function(df, grp) {
    lin_mod  <- fit_linear(df)
    exp_mod  <- fit_exponen(df)
    shift    <- min(df$catch_monthly_mean) - 0.1
    mnths    <- seq(1, 12, by = 0.2)
    lin_pred <- predict(lin_mod, newdata = data.frame(month_num = mnths))
    exp_pred <- if (!is.null(exp_mod))
      exp(predict(exp_mod, newdata = data.frame(month_num = mnths))) + shift
    else NA_real_
    tibble(year = grp$year, month_num = mnths,
           linear = lin_pred, exponential = exp_pred)
  }) %>%
  bind_rows() %>%
  pivot_longer(c(linear, exponential), names_to = "model", values_to = "pred_temp")

p3 <- ggplot() +
  geom_ribbon(data = catch_plot,
              aes(x = month_num,
                  ymin = catch_monthly_mean - catch_monthly_sd,
                  ymax = catch_monthly_mean + catch_monthly_sd,
                  fill = year_f), alpha = 0.15) +
  geom_line(data = catch_plot,
            aes(x = month_num, y = catch_monthly_mean,
                colour = year_f), linewidth = 1.2) +
  geom_point(data = catch_plot,
             aes(x = month_num, y = catch_monthly_mean,
                 colour = year_f), size = 2.5) +
  geom_line(data = pred_data,
            aes(x = month_num, y = pred_temp,
                linetype = model, colour = factor(year)),
            linewidth = 0.9, alpha = 0.85) +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb) +
  scale_colour_manual(values = c("2023" = "#E07B39", "2024" = "#3B7EC8"),
                      name = "Year") +
  scale_fill_manual(values   = c("2023" = "#E07B39", "2024" = "#3B7EC8"),
                    name = "Year (±1 SD)") +
  scale_linetype_manual(values = c("linear" = "dashed", "exponential" = "dotdash"),
                        name = "Trend model") +
  labs(title   = "Catchment Monthly Mean Temperature with Trend Models",
       subtitle = "Solid line = observed | Dashed = linear fit | Dot-dash = exponential fit\nRibbon = ± 1 SD across daily observations",
       x = "Month", y = "Mean Temperature (°C)") +
  theme_boye

ggsave(file.path(plots_dir, "03_catchment_trend.png"), p3,
       width = 13, height = 7, dpi = 180)
cat("Saved plot 3\n")

# ── 7.4 Daily temperature range boxplots by month ─────────────────────────────
p4 <- daily_stats %>%
  ggplot(aes(x = month_label, y = daily_range, fill = factor(year))) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.75, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("2023" = "#E07B39", "2024" = "#3B7EC8"), name = "Year") +
  labs(title   = "Daily Temperature Range Distribution by Month",
       subtitle = "Box = interquartile range | Line = median | Whiskers = 1.5 × IQR",
       x = NULL, y = "Daily Temperature Range (°C)") +
  theme_boye

ggsave(file.path(plots_dir, "04_daily_range_boxplots.png"), p4,
       width = 13, height = 6, dpi = 180)
cat("Saved plot 4\n")

# ── 7.5 Yearly comparison: 2023 vs 2024 mean temp per stream ──────────────────
p5 <- yearly_stats %>%
  ggplot(aes(x = reorder(stream, yearly_mean), y = yearly_mean,
             fill = factor(year))) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.9) +
  geom_errorbar(aes(ymin = yearly_mean - yearly_sd,
                    ymax = yearly_mean + yearly_sd),
                position = position_dodge(0.75), width = 0.3, linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("2023" = "#E07B39", "2024" = "#3B7EC8"), name = "Year") +
  labs(title   = "Annual Mean Temperature per Stream – 2023 vs 2024",
       subtitle = "Error bars = ± 1 SD across monthly means",
       x = NULL, y = "Annual Mean Temperature (°C)") +
  theme_boye

ggsave(file.path(plots_dir, "05_yearly_comparison.png"), p5,
       width = 11, height = 9, dpi = 180)
cat("Saved plot 5\n")

# ── 7.6 Per-stream monthly trends with linear & exponential fits ──────────────
# Build prediction ribbons for each stream
stream_pred <- monthly_stats %>%
  group_by(stream, year) %>%
  group_map(function(df, grp) {
    if (nrow(df) < 3) return(NULL)
    shift <- min(df$monthly_mean) - 0.1
    df2   <- df %>% mutate(temp_pos = monthly_mean - shift,
                           month_num = as.numeric(month))
    lin_m <- tryCatch(lm(monthly_mean ~ month_num, data = df2), error = function(e) NULL)
    exp_m <- tryCatch(lm(log(temp_pos) ~ month_num, data = df2), error = function(e) NULL)
    mnths <- seq(min(df$month), max(df$month), by = 0.3)
    nd    <- data.frame(month_num = mnths)
    lin_p <- if (!is.null(lin_m)) predict(lin_m, nd) else NA_real_
    exp_p <- if (!is.null(exp_m)) exp(predict(exp_m, nd)) + shift else NA_real_
    tibble(stream = grp$stream, year = grp$year, month_num = mnths,
           linear = lin_p, exponential = exp_p)
  }) %>%
  bind_rows() %>%
  pivot_longer(c(linear, exponential), names_to = "model", values_to = "pred")

p6 <- ggplot() +
  geom_line(data = monthly_stats,
            aes(x = month, y = monthly_mean, colour = factor(year)), linewidth = 0.8) +
  geom_point(data = monthly_stats,
             aes(x = month, y = monthly_mean, colour = factor(year)), size = 1.5) +
  geom_line(data = stream_pred,
            aes(x = month_num, y = pred,
                colour = factor(year), linetype = model), linewidth = 0.6, alpha = 0.8) +
  facet_wrap(~ stream, ncol = 4, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan","Apr","Jul","Oct")) +
  scale_colour_manual(values = c("2023" = "#E07B39", "2024" = "#3B7EC8"), name = "Year") +
  scale_linetype_manual(values = c("linear" = "dashed", "exponential" = "dotdash"),
                        name = "Trend") +
  labs(title   = "Monthly Temperature Trends per Stream – Linear & Exponential Fits",
       subtitle = "Solid = observed | Dashed = linear | Dot-dash = exponential",
       x = "Month", y = "Mean Temperature (°C)") +
  theme_boye +
  theme(axis.text.x = element_text(size = 7))

ggsave(file.path(plots_dir, "06_stream_trends.png"), p6,
       width = 18, height = 14, dpi = 180)
cat("Saved plot 6\n")

# ── 7.7 Catchment daily mean temperature ribbon (annual cycle) ─────────────────
p7 <- catch_daily %>%
  mutate(doy = yday(date)) %>%
  group_by(year, doy) %>%
  summarise(mean_t = mean(catch_mean, na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = doy, y = mean_t, colour = factor(year), fill = factor(year))) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, alpha = 0.15) +
  scale_x_continuous(
    breaks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)),
    labels = month.abb) +
  scale_colour_manual(values = c("2023" = "#E07B39", "2024" = "#3B7EC8"), name = "Year") +
  scale_fill_manual(values   = c("2023" = "#E07B39", "2024" = "#3B7EC8"), name = "Year") +
  labs(title   = "Catchment Daily Mean Temperature – Annual Cycle",
       subtitle = "Shaded band = LOESS smooth ± 95% CI",
       x = "Day of Year (month labels at month start)",
       y = "Catchment Mean Temperature (°C)") +
  theme_boye

ggsave(file.path(plots_dir, "07_catchment_annual_cycle.png"), p7,
       width = 14, height = 6, dpi = 180)
cat("Saved plot 7\n")

cat("\nAll 7 plots saved to:", plots_dir, "\n")

# =============================================================================
# SECTION 8 – Print summaries to console
# =============================================================================
cat("\n══════════════════════════════════════════════\n")
cat("  CATCHMENT MONTHLY SUMMARY\n")
cat("══════════════════════════════════════════════\n")
print(catch_monthly %>% select(-catch_monthly_sd) %>% as.data.frame(), row.names = FALSE)

cat("\n══════════════════════════════════════════════\n")
cat("  CATCHMENT YEARLY SUMMARY\n")
cat("══════════════════════════════════════════════\n")
print(catch_yearly, row.names = FALSE)

cat("\n══════════════════════════════════════════════\n")
cat("  YEARLY STATS PER STREAM\n")
cat("══════════════════════════════════════════════\n")
print(yearly_stats %>% as.data.frame(), row.names = FALSE)

# =============================================================================
# SECTION 9 – OPTIONAL: Save all results to a single Excel workbook
# =============================================================================
#Uncomment the block below to export all tables to one Excel file.
# Requires the 'openxlsx' package:
#   install.packages("openxlsx")
#
# ─────────────────────────────────────────────────────────────────────────────
library(openxlsx)

out_excel <- file.path(folder, "Boye_Temperature_Analysis_2023_2024.xlsx")

wb <- createWorkbook()

hs <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#1F497D",
  halign = "CENTER",
  textDecoration = "Bold"
)

sheets_list <- list(
  metadata          = metadata,
  daily_stats       = daily_stats,
  monthly_stats     = monthly_stats,
  yearly_stats      = yearly_stats,
  catchment_daily   = catch_daily,
  catchment_monthly = catch_monthly,
  catchment_yearly  = catch_yearly
)

for (nm in names(sheets_list)) {
  
  addWorksheet(wb, nm)
  
  writeData(
    wb,
    nm,
    sheets_list[[nm]],
    headerStyle = hs
  )
  
  setColWidths(
    wb,
    nm,
    cols = seq_len(ncol(sheets_list[[nm]])),
    widths = "auto"
  )
  
  addFilter(
    wb,
    nm,
    row = 1,
    cols = seq_len(ncol(sheets_list[[nm]]))
  )
  
}

for (nm in names(sheets_list)) {
  freezePane(wb, nm, firstRow = TRUE)
}

saveWorkbook(wb, out_excel, overwrite = TRUE)

cat("\nExcel workbook saved to:", out_excel, "\n")
# ─────────────────────────────────────────────────────────────────────────────

cat("\nAnalysis complete.\n")
cat("Data frames in memory: raw_all, daily_stats, monthly_stats, yearly_stats,\n")
cat("                        catch_daily, catch_monthly, catch_yearly, metadata\n")
