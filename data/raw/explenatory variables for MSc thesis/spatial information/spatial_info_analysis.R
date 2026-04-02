library(dplyr)

# ── Name lookup tables ────────────────────────────────────────────────────────

boye_name_lookup <- c(
  "BOY_kl"       = "BOYkl",
  "BOY_oh_B224"  = "BOYohB224",
  "BOY_oh_Br"    = "BOYohBr",
  "BOY_oh_Ki"    = "BOYohKi",
  "BOY_oh_Sp"    = "BOYohSp",
  "BOY_uh_Ha"    = "BOYuhHa",
  "BOY_uh_Sp"    = "BOYuhSp",
  "BRA_ob"       = "BRAob",
  "BRA_oh_Bo"    = "BRAohBo",
  "HAA_ob"       = "HAAob",
  "HAA_un"       = "HAAun",
  "KIR_ob"       = "KIRob",
  "KIR_un"       = "KIRun",
  "NAT_oh_Bo"    = "NATohBo",
  "Qu_oh_Bo"     = "QUAohBo",
  "SCH_oh_Vo"    = "SCHohVo",
  "VOR_oh_Bo"    = "VORohBo",
  "VOR_uh_Sch"   = "VORuhSc",
  "WIT_ob"       = "WITob"
)

rotbach_name_lookup <- c(
  "Rotbach"     = "ROTB",
  "Schwarzbach" = "SCHW"
)

# ── Load and clean each file ──────────────────────────────────────────────────

boye <- read.csv("Boye_all_sites_croppercent.csv") %>%
  mutate(
    stream = "Boye",
    site   = recode(Name, !!!boye_name_lookup)
  )

kinzig <- read.csv("Kinzig_all_sites_croppercent.csv") %>%
  mutate(
    stream = "Kinzig",
    site   = Name
  )

rotbach <- read.csv("Rotbach_Schwarzbach_Landuse.csv") %>%
  mutate(
    stream          = "Boye",
    site            = recode(Name, !!!rotbach_name_lookup),
    Area_Percentage = as.numeric(gsub("%", "", Area_Percentage))  # fix character type
  )

# ── Combine ───────────────────────────────────────────────────────────────────

all_sites <- bind_rows(boye, kinzig, rotbach) %>%
  select(stream, site, Year, Value, Area, Area_Percentage)







# Quick overview
glimpse(all_sites)

# Check all unique site names are clean
all_sites %>% distinct(stream, site) %>% arrange(stream, site) %>% print(n = Inf)

# Confirm no NAs introduced in Area_Percentage (from the % stripping)
all_sites %>% filter(is.na(Area_Percentage))

# Row counts per source — sanity check
all_sites %>% count(stream)

# First few rows
head(all_sites)


# Option 1 — convert to tibble first, then print all rows
all_sites %>% distinct(stream, site) %>% arrange(stream, site) %>% as_tibble() %>% print(n = Inf)

# Option 2 — base R, always works regardless of class
all_sites %>% distinct(stream, site) %>% arrange(stream, site) %>% as.data.frame()

# Option 3 — open in RStudio's viewer pane (most convenient)
all_sites %>% distinct(stream, site) %>% arrange(stream, site) %>% View()

View(all_sites)

library(dplyr)
library(openxlsx)  # install.packages("openxlsx") if needed

# ── Filter merged data to 2022 only ──────────────────────────────────────────

data_2022 <- all_sites %>%
  filter(Year == 2022)

# ── Create workbook with two sheets ──────────────────────────────────────────

wb <- createWorkbook()

# Sheet 1: full 2022 data
addWorksheet(wb, "All_Sites_2022")
writeDataTable(wb, sheet = "All_Sites_2022", x = data_2022, tableStyle = "TableStyleMedium9")

# Sheet 2: summary — top land use per site
summary_2022 <- data_2022 %>%
  group_by(stream, site) %>%
  slice_max(Area_Percentage, n = 1) %>%
  ungroup() %>%
  arrange(stream, site)

addWorksheet(wb, "Summary_2022")
writeDataTable(wb, sheet = "Summary_2022", x = summary_2022, tableStyle = "TableStyleMedium2")

# ── Style: widen columns for readability ─────────────────────────────────────

setColWidths(wb, "All_Sites_2022", cols = 1:ncol(data_2022),    widths = "auto")
setColWidths(wb, "Summary_2022",   cols = 1:ncol(summary_2022), widths = "auto")

# ── Save ──────────────────────────────────────────────────────────────────────

saveWorkbook(wb, "Landuse_2022_merged.xlsx", overwrite = TRUE)

message("Saved: Landuse_2022_merged.xlsx")


library(dplyr)

# ── Land use category definitions ────────────────────────────────────────────

agricultural_codes <- c(
  115,  # Winter wheat
  121,  # Winter rye
  131,  # Winter barley
  132,  # Spring barley
  140,  # Oats
  155,  # Triticale
  171,  # Grain corn
  311,  # Winter rape
  411,  # Silage corn
  422,  # Clover grass
  424,  # Grass
  459,  # Grass (pasture)
  602,  # Potatoes
  9012, # Non-irrigated arable land
  9018, # Pastures
  9020, # Complex cultivation patterns
  9021  # Land principally occupied by agriculture
)

industrial_codes <- c(
  9003, # Industrial or commercial units
  9007, # Mineral extraction sites
  9008  # Dump sites
)

urban_codes <- c(
  9001, # Continuous urban fabric
  9002, # Discontinuous urban fabric
  9010, # Green urban areas
  9011  # Sport and leisure facilities
)

rural_codes <- c(
  9023, # Broad-leaved forest
  9024, # Coniferous forest
  9025, # Mixed forest
  9029, # Transitional woodland-shrub
  9040, # Water courses
  9041  # Water bodies
)

# ── Tag each row with its category ───────────────────────────────────────────

data_2022_cat <- data_2022 %>%
  mutate(
    category = case_when(
      Value %in% agricultural_codes ~ "agricultural",
      Value %in% industrial_codes   ~ "industrial",
      Value %in% urban_codes        ~ "urban",
      Value %in% rural_codes        ~ "rural",
      TRUE                          ~ "other"
    )
  )

# ── Build summary table ───────────────────────────────────────────────────────

summary_2022 <- data_2022_cat %>%
  group_by(stream, site) %>%
  summarise(
    total_area             = sum(Area, na.rm = TRUE),
    total_area_agric       = sum(Area[category == "agricultural"], na.rm = TRUE),
    total_area_industrial  = sum(Area[category == "industrial"],   na.rm = TRUE),
    total_area_urban       = sum(Area[category == "urban"],        na.rm = TRUE),
    total_area_rural       = sum(Area[category == "rural"],        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_agricultural  = round(total_area_agric      / total_area * 100, 2),
    pct_industrial    = round(total_area_industrial  / total_area * 100, 2),
    pct_urban         = round(total_area_urban       / total_area * 100, 2),
    pct_rural         = round(total_area_rural       / total_area * 100, 2)
  ) %>%
  arrange(stream, site)

# ── Quick check: do categories cover all area? ───────────────────────────────

summary_2022 %>%
  mutate(pct_total = pct_agricultural + pct_industrial + pct_urban + pct_rural) %>%
  select(stream, site, pct_total) %>%
  filter(pct_total < 95)   # flag sites where >5% fell into "other"

# ── View result ───────────────────────────────────────────────────────────────

View(summary_2022)

# ── Check: flag anything that fell into "other" ───────────────────────────────

data_2022_cat %>%
  filter(category == "other") %>%
  distinct(Value)   # should return empty

# ── Save to Excel ─────────────────────────────────────────────────────────────

wb <- createWorkbook()
addWorksheet(wb, "Summary_2022")
writeDataTable(wb, sheet = "Summary_2022", x = summary_2022, tableStyle = "TableStyleMedium9")
setColWidths(wb, "Summary_2022", cols = 1:ncol(summary_2022), widths = "auto")

saveWorkbook(wb, "LandUse_Summary_2022.xlsx", overwrite = TRUE)
message("Saved: LandUse_Summary_2022.xlsx")