# ===============================
# FIELD DATE SUMMARY TO EXCEL
# ===============================

# Run this ONCE to install packages (remove # to uncomment):
install.packages(c("readr", "dplyr", "lubridate", "openxlsx"))

library(readr)
library(dplyr)
library(lubridate)
library(openxlsx)

# 1. Read the file
df <- read_csv("2025_field_processed_samples - Copy (2).csv")

# 2. Convert dates
df <- df %>%
  mutate(
    date_deployment = dmy(date_deployment),
    date_sampled    = dmy(date_sampled)
  )

# 3. Keep only Boye and Kinzig
df_bk <- df %>%
  filter(catchment %in% c("Boye", "Kinzig"))

# 4. Add exact incubation days per bag
df_bk <- df_bk %>%
  mutate(
    incubation_days = as.numeric(date_sampled - date_deployment)
  )

# 5. Create catchment summary table
summary_table <- df_bk %>%
  group_by(catchment) %>%
  summarise(
    first_deployment       = min(date_deployment, na.rm = TRUE),
    last_deployment        = max(date_deployment, na.rm = TRUE),
    first_sample           = min(date_sampled, na.rm = TRUE),
    last_sample            = max(date_sampled, na.rm = TRUE),
    deployment_window_days = as.numeric(last_deployment - first_deployment),
    sampling_window_days   = as.numeric(last_sample - first_sample),
    min_incubation_days    = min(incubation_days, na.rm = TRUE),
    max_incubation_days    = max(incubation_days, na.rm = TRUE)
  ) %>%
  ungroup()

# 6. Create a nice text summary table
# FIX: deployment_range and sampling_range are now included in select()
summary_text <- summary_table %>%
  mutate(
    deployment_range = paste0(format(first_deployment, "%d %b %Y"),
                              " to ",
                              format(last_deployment, "%d %b %Y")),
    sampling_range   = paste0(format(first_sample, "%d %b %Y"),
                              " to ",
                              format(last_sample, "%d %b %Y")),
    incubation_range = paste0(min_incubation_days, "–", max_incubation_days, " days")
  ) %>%
  select(
    catchment,
    deployment_range,      # ← was missing
    sampling_range,        # ← was missing
    deployment_window_days,
    sampling_window_days,
    incubation_range
  )

# 7. Detailed per-bag table
bag_details <- df_bk %>%
  select(
    catchment, stream, type, replicate, bag_ID,
    date_deployment, date_sampled, incubation_days
  ) %>%
  arrange(catchment, stream, type, replicate)

# 8. Catchment-specific tables
boye_details   <- bag_details %>% filter(catchment == "Boye")
kinzig_details <- bag_details %>% filter(catchment == "Kinzig")

# 9. Write to Excel
wb <- createWorkbook()

addWorksheet(wb, "Summary");        writeData(wb, "Summary",        summary_table)
addWorksheet(wb, "Summary_Text");   writeData(wb, "Summary_Text",   summary_text)
addWorksheet(wb, "Bag_Details_All"); writeData(wb, "Bag_Details_All", bag_details)
addWorksheet(wb, "Boye_Details");   writeData(wb, "Boye_Details",   boye_details)
addWorksheet(wb, "Kinzig_Details"); writeData(wb, "Kinzig_Details", kinzig_details)

for (sheet in names(wb)) {
  setColWidths(wb, sheet = sheet, cols = 1:20, widths = "auto")
}

saveWorkbook(wb, "Boye_Kinzig_Field_Date_Summary.xlsx", overwrite = TRUE)

# 10. Print result
print(summary_table)
print(summary_text)
cat("\nExcel file saved as: Boye_Kinzig_Field_Date_Summary.xlsx\n")