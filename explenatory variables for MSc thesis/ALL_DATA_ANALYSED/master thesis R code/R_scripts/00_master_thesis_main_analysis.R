# =============================================================================
# The Effect of Multiple Stressors on Microbial and Total Leaf Decomposition
# in Urban and Rural Streams of the Kinzig Catchment and
# Restored and Unrestored Streams of the Boye Catchment
#
# Author:     Ikenna Stephen Agube
# Supervisor: Dr. Verena C. Schreiner
# Institute:  Department of Ecotoxicology, University of Duisburg-Essen
# Contact:    ikenna.agube@uni-due.de;
# =============================================================================
#
# NOTE ON EXPERIMENTAL DESIGN:
#
# Kinzig catchment – predictor: land use (rural vs. urban)
#   Research question: Does urbanisation reduce leaf litter decomposition?
# Boye catchment – predictor: restoration status
#   (reference / recently restored ~1 yr / long restored ~20 yr)
#   Research question: Does stream restoration recover leaf litter decomposition,
#   and does recovery time matter for microbial vs. macroinvertebrate-mediated
#   breakdown?
#
# Both catchments share the same raw data files (processed, initial_mass,
# transport, inverts); catchment-specific physico-chemistry is in separate files.
# =============================================================================


# =============================================================================
# STUDY HYPOTHESES
# =============================================================================
#
# H1 – GENERAL HYPOTHESIS (both catchments)
#   A higher proportion of agricultural land use in the upstream catchment, as
#   well as elevated conductivity, negatively affects leaf decomposition across
#   both the Boye and Kinzig catchments. Agricultural land has been identified
#   as a key negative predictor of microbial decomposition (Schreiner et al.,
#   2023), and conductivity reflects the degree of anthropogenic alteration
#   associated with agricultural and urban-industrial land use (Bundschuh et
#   al., 2021).
#   --> Tested in Section 8: lmer(k_microbial ~ pct_agri + conductivity +
#       (1|stream), data = all_sites)
#
# H2 – COMPARATIVE HYPOTHESIS (Boye vs Kinzig)
#   Microbial leaf decomposition is differently affected by environmental
#   stressors in the two catchments. In the Boye, higher toxic units of trace
#   metals (Cd, Cu, Fe, Ni, Pb, Zn) are associated with reduced microbial
#   decomposition, reflecting the sensitivity of aquatic hyphomycetes and
#   bacteria to metal contamination in urban-industrial streams (Duarte et al.,
#   2008; Ferreira et al., 2016; Feyen et al., 2022). In the Kinzig, higher
#   nitrate and phosphate concentrations are associated with enhanced microbial
#   decomposition, as agricultural drainage can enrich streams with nutrients
#   that stimulate fungal and bacterial decomposer activity (Gulis &
#   Suberkropp, 2003; Robinson & Gessner, 2000).
#   --> Tested in Section 9: separate lmer models per catchment, with
#       metal_TU (Boye) and nutrient concentrations (Kinzig) as predictors
#
# H3 – FINE SUBSTRATE (Boye)
#   Sites in the Boye catchment with higher proportions of fine substrate
#   (FPOM + sand, >= 25% streambed cover) exhibit slower total leaf
#   decomposition (k_total). Macroinvertebrate shredders, which typically
#   drive 60-80% of mass loss in coarse-substrate streams, become physically
#   excluded when fine sediments bury leaf packs and eliminate accessible
#   habitat (Lepori et al., 2005). Microbial decomposers may benefit from
#   prolonged leaf-water contact (Pascoal & Cassio, 2004), but microbes alone
#   decompose at 0.002-0.006 day-1 (Suter et al., 2011), which is
#   insufficient to offset shredder loss.
#   --> Tested in Section 10: high_fine (>= 25%) vs low_fine (< 25%) groups,
#       lmer(k_total ~ fine_class + (1|stream), data = boye_total)
#
# H4 – COARSE SUBSTRATE (Kinzig)
#   Greater availability of coarse particulate organic matter (CPOM) in the
#   Kinzig catchment positively enhances total leaf decomposition. Higher CPOM
#   availability provides more habitat and food for shredding macroinvertebrates,
#   supporting a more diverse and abundant shredder community and thereby
#   increasing the macroinvertebrate contribution to overall leaf mass loss
#   (Graca et al., 2001).
#   --> Tested in Section 11: lmer(k_total ~ cpom_pct + (1|stream),
#       data = kinzig_total)
#
# =============================================================================


# 1 Load tools / packages ------------------------------------------------------

library(tidyverse)      # version 2.0.0, data manipulation and plotting
library(readxl)         # version 1.4.3, read Excel files
library(lme4)           # version 1.1-35, fit mixed-effects models
library(lmerTest)       # version 3.1-3, p-values for lme4 models
library(emmeans)        # version 1.10, estimated marginal means + pairwise contrasts
library(performance)    # version 0.11, model diagnostics
library(janitor)        # version 2.2.0, clean column names
library(ggeffects)      # version 1.3.4, partial effects plots
library(cowplot)        # version 1.1.3, arrange multi-panel figures
library(broom.mixed)    # version 0.2.9, tidy mixed model output
library(stringr)        # version 1.5.1, string manipulation
library(glue)           # version 1.7.0, formatted messages
library(e1071)          # skewness() for distribution diagnostics (Schreiner Step 2)
library(car)            # vif() for collinearity diagnostics  (Schreiner Step 3)
library(glmnet)         # LASSO fitting, required by stabs    (Schreiner Step 4)
library(stabs)          # stabsel() stability selection        (Schreiner Step 4)
library(ggsignif)       # geom_signif() significance brackets on plots
library(multcomp)       # cld() compact letter display for post-hoc groups
library(multcompView)   # required by emmeans 2.x for cld() letter display

options(scipen = 999, digits = 3)


# set path – points to data/processed/ on your Mac
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/masters/Master thesis documents to use/explenatory variables for MSc thesis/ALL_DATA_ANALYSED/master thesis R code/data/processed")


# 2 Import data ----------------------------------------------------------------

# All files live in data/processed/ relative to the project root.
# The working directory is set above to that folder, so file names alone are
# sufficient for read_excel() and read_csv().

## 2.1 Leaf-bag field data (shared across both catchments) ---------------------

# Main processed leaf-bag data: one row per bag retrieved from the field.
# Used in sections 3–9 for AFDM calculation, k computation, and models.
processed <- read_excel("2025_field_processed_samples.xlsx")

# Transport controls: bags that never entered the stream, used to quantify
# handling loss (section 3.3). Columns: bag_ID, type, catchment,
# initial_mass_g, mass_60_g, mass_500_g, comment.
transport <- read_excel("2025_field_transport_controls.xlsx")

# Starting dry masses for every bag before deployment (section 3.4).
# Columns: bag_ID, type, initial_mass_g, comment.
initial_mass <- read_excel("2025_field_initial_mass.xlsx")


## 2.2 Kinzig-specific files ---------------------------------------------------

# Site-level land use classification for the Kinzig catchment (section 4.1).
# Columns: stream, catchment, position, rural, urban.
site_info_kinzig <- read_excel("2025_field_site_information_Kinzig.xlsx")

# Physico-chemical measurements for Kinzig sites (section 4.2).
# Raw column names contain newlines / special characters; clean_names() is
# applied in section 4.2, which produces a 'water_temp' column.
phch_kinzig <- read_excel("202505_phch_data_Kinzig.xlsx")


## 2.3 Boye-specific files -----------------------------------------------------

# Water-chemistry data for Boye sites (section 5).
# clean_names() converts 'water_temp_degC' → 'water_temp_deg_c'.
# Section 5.2 then renames site → stream and water_temp_deg_c → water_temp.
phch_boye <- read_excel("boye_processed_samples_mar2025.xlsx") %>%
  clean_names() %>%
  # Many measurement columns arrive as character (e.g. "<LOD" entries in the
  # Excel sheet).  Coerce all columns except the site identifier to numeric;
  # values that cannot be parsed (like "<LOD") become NA silently.
  mutate(across(-site, ~suppressWarnings(as.numeric(.))))

# Temperature logger: one row per date at each Boye site.
# Used in section 5.2 to calculate mean water temperature and degree days.
temp_logger_boye <- read_csv("boye_temperature_logger_march2025.csv",
                             show_col_types = FALSE)


## 2.4 Master environmental variables (both catchments) ------------------------

# Comprehensive site-level table covering all 44 streams in Boye and Kinzig.
# Contains physicochemistry (pH, conductivity, DO, nutrients, metals),
# mean water temperature from loggers (Avg_Temp_degC), substrate composition
# (Makro_pct … Tech_pct), riparian vegetation (Submers_pct … FPOM_pct),
# and catchment land-use areas and percentages (Pct_Agri, Pct_Urban, etc.).
# clean_names() standardises all column headers to snake_case.
#
# Columns after clean_names():
#   stream, catchment,
#   ph, conductivity_u_s_cm, do_mg_l, o2_sat_pct, water_temp_deg_c,
#   water_level_cm, ortho_po4_mg_l, total_po4_mg_l, nh4_mg_l, no2_mg_l,
#   no3_mg_l, tn_mg_l, cl_mg_l, so4_mg_l, avg_temp_deg_c,
#   cd_mg_l, cu_mg_l, fe_mg_l, ni_mg_l, pb_mg_l, zn_mg_l,
#   ca_mg_l, mg_mg_l, na_mg_l, k_mg_l, co3_mg_l, hco3_mg_l,
#   makro_pct, meso_pct, mikro_pct, akal_pct, psam_pct, tech_pct,
#   submers_pct, emers_pct, living_parts_pct, xylal_pct, cpom_pct, fpom_pct,
#   total_area_km2, area_agric, area_indus, area_urban, area_rural,
#   pct_agri, pct_indus, pct_urban, pct_rural, comments, date

master_vars <- read_excel("Master_DAta_VAriables_RESIST_Field_MArch_2025.xlsx") %>%
  clean_names() %>%
  mutate(catchment = str_to_title(str_trim(catchment)),   # "Boye" / "Kinzig"
         stream    = str_trim(stream)) %>%
  # Many numeric columns arrive as character because Excel cells contain text
  # like "n.a.", "-", or "<LOD".  Coerce all character columns except the
  # string identifiers and free-text fields to numeric; non-parseable entries
  # silently become NA.
  mutate(across(where(is.character) & !all_of(c("catchment", "stream", "comments")),
                ~suppressWarnings(as.numeric(.))))

# Split into catchment-specific subsets for convenient use downstream
master_boye   <- master_vars %>% filter(catchment == "Boye")
master_kinzig <- master_vars %>% filter(catchment == "Kinzig")

message(glue(
  "master_vars loaded: {nrow(master_vars)} sites  ",
  "({nrow(master_boye)} Boye, {nrow(master_kinzig)} Kinzig)"
))


# 3 Shared pre-processing: AFDM correction for all catchments ------------------

## 3.1 Exposure duration -------------------------------------------------------

# Convert deployment and sampling dates to Date objects and compute deployment
# duration in days. This is used later when calculating degree days.

processed$date_deployment <- as.Date(processed$date_deployment, format = "%d/%m/%Y")
processed$date_sampled    <- as.Date(processed$date_sampled,    format = "%d/%m/%Y")
processed$deployment_duration_days <- as.numeric(
  difftime(processed$date_sampled, processed$date_deployment, units = "days")
)


## 3.2 Clean variable types ----------------------------------------------------

# Ensure numeric columns are truly numeric. Files with missing values often
# import numeric fields as characters; this loop coerces them safely.

numeric_cols <- c("OMD_60_mg", "OMD_500_mg", "no_disks_OMD", "no_disks_all",
                  "rest_60_g", "rest_500_g", "initial_mass_g")
for (col in numeric_cols) {
  if (col %in% names(processed)) {
    processed[[col]] <- as.numeric(processed[[col]])
  }
}

# Standardise type labels (fine / coarse) to lower case
processed <- processed %>% mutate(type = tolower(str_trim(type)))


## 3.3 Handling-loss correction ------------------------------------------------

# Ash-free dry mass (AFDM) corrects for inorganic material on the leaf bags.
# We calculate the AFDM fraction lost between 60°C and 500°C in the transport
# controls, separately for fine and coarse mesh types.

transport <- transport %>%
  mutate(catchment = tolower(str_trim(catchment)),
         type      = tolower(str_trim(type))) %>%
  mutate(
    perc_60   = mass_60_g  / initial_mass_g,
    perc_500  = mass_500_g / initial_mass_g,
    afdm_frac = perc_60 - perc_500
  ) %>%
  filter(is.finite(afdm_frac), afdm_frac > 0, afdm_frac < 1)

# Mean AFDM fraction per catchment and mesh type
corr_afdm <- transport %>%
  group_by(catchment, type) %>%
  summarise(afdm_frac = mean(afdm_frac, na.rm = TRUE), .groups = "drop")

stopifnot(nrow(corr_afdm) >= 1)


## 3.4 Join initial masses and compute AFDM start ------------------------------

# Some rows may lack initial_mass_g because the values reside in a separate
# file. We match by bag_ID and type.

if (!"initial_mass_g" %in% names(processed)) {
  initial_mass <- initial_mass %>% mutate(type = tolower(str_trim(type)))
  processed <- processed %>%
    left_join(initial_mass %>% select(bag_ID, type, initial_mass_g),
              by = c("bag_ID", "type"))
}

# Join handling-loss fractions. The catchment column in processed must match
# the catchment column in corr_afdm (both lower case).
processed <- processed %>%
  mutate(catchment_lc = tolower(str_trim(catchment))) %>%
  left_join(corr_afdm, by = c("catchment_lc" = "catchment", "type")) %>%
  mutate(AFDM_start = initial_mass_g * afdm_frac) %>%
  relocate(AFDM_start, .after = initial_mass_g)


## 3.5 Disk correction and AFDM end --------------------------------------------

# In a companion study, leaf disks were cut from the bags for additional
# endpoints. The AFDM of these disks was also determined; the remaining
# organic matter is now corrected for these removed disks. Where bag-level
# information is unavailable, we fall back to the stream mean and then to
# the global mean. Masses are converted from mg to g.

processed <- processed %>%
  mutate(no_disks_OMD = replace_na(no_disks_OMD, 0),
         no_disks_all = replace_na(no_disks_all, 0)) %>%
  mutate(
    mean_afdm_per_disk = ifelse(
      is.finite(OMD_60_mg) & is.finite(OMD_500_mg) & no_disks_OMD > 0,
      (OMD_60_mg - OMD_500_mg) / no_disks_OMD,
      NA_real_)
  )

# stream-level mean AFDM per disk
stream_means <- processed %>%
  group_by(stream) %>%
  summarise(stream_mean_disk = mean(mean_afdm_per_disk, na.rm = TRUE),
            .groups = "drop")

# global fallback (in mg)
global_mean <- processed %>%
  filter(is.finite(mean_afdm_per_disk)) %>%
  summarise(val = mean(mean_afdm_per_disk, na.rm = TRUE)) %>%
  pull(val)
if (!is.finite(global_mean)) global_mean <- 0

# total AFDM end (disk AFDM converted from mg to g)
processed <- processed %>%
  left_join(stream_means, by = "stream") %>%
  mutate(
    mean_afdm_per_disk = coalesce(mean_afdm_per_disk, stream_mean_disk, global_mean),
    afdm_all_disks     = (no_disks_all * mean_afdm_per_disk) / 1000,
    rest_AFDM          = rest_60_g - rest_500_g,
    AFDM_end           = rest_AFDM + afdm_all_disks
  ) %>%
  select(-stream_mean_disk)

# Remove rows where core mass values are missing
n0 <- nrow(processed)
processed <- processed %>%
  filter(is.finite(rest_60_g), is.finite(rest_500_g),
         is.finite(AFDM_start), is.finite(AFDM_end))
message(glue("Dropped {n0 - nrow(processed)} rows due to missing masses."))


# 4 Kinzig catchment -----------------------------------------------------------

## 4.1 Land use assignment -----------------------------------------------------

# Translate the land use flags from site_info into a factor with levels
# "rural" and "urban". Rural is set as the reference level.

site_info_kinzig <- site_info_kinzig %>%
  mutate(LandUse = case_when(
    tolower(str_trim(rural)) == "x" ~ "rural",
    tolower(str_trim(urban)) == "x" ~ "urban",
    TRUE ~ NA_character_
  )) %>%
  select(stream, LandUse)


## 4.2 Degree days and k – Kinzig ----------------------------------------------

# Summarise Kinzig physico-chemistry by stream to get mean water temperature.

phch_kinzig_summ <- phch_kinzig %>%
  clean_names() %>%
  group_by(stream) %>%
  summarise(water_temp = mean(water_temp, na.rm = TRUE), .groups = "drop")

# Filter to Kinzig and compute k per degree day:
#
#   k = -ln(AFDM_t1 / AFDM_t0) / (T * n)
#
# where AFDM_t0 is the initial AFDM, AFDM_t1 the AFDM after deployment,
# T is the mean stream temperature and n the number of deployment days.

kinzig <- processed %>%
  filter(catchment == "Kinzig") %>%
  left_join(site_info_kinzig, by = "stream") %>%
  left_join(phch_kinzig_summ,  by = "stream") %>%
  mutate(
    deg_days = deployment_duration_days * water_temp,
    k        = -log(AFDM_end / AFDM_start) / deg_days
  ) %>%
  filter(is.finite(k), k > 0, AFDM_end < AFDM_start, deg_days > 0,
         !is.na(LandUse)) %>%
  mutate(
    k       = round(k, 5),
    LandUse = factor(LandUse, levels = c("rural", "urban"))
  )

message(glue("Kinzig: {nrow(kinzig)} bags after filtering."))


## 4.3 Join master environmental predictors – Kinzig ---------------------------
#
# We bring in the key predictor columns from master_kinzig for each stream:
# agricultural land use percentage, conductivity, NO3, ortho-PO4, and CPOM
# coverage. These are used in the hypothesis-specific models (H1, H2, H4).
# We use a left_join so that bags without a matching master row are kept but
# will have NA predictors, which the models will silently drop.

kinzig <- kinzig %>%
  left_join(
    master_kinzig %>%
      select(stream,
             pct_agri,
             conductivity_u_s_cm,
             no3_mg_l,
             ortho_po4_mg_l,
             total_po4_mg_l,
             cpom_pct,
             fpom_pct,
             psam_pct),
    by = "stream"
  )


## 4.4 Prepare Kinzig analysis datasets -----------------------------------------
#
# Separate microbial (fine mesh) and total (coarse mesh) datasets.
# These are used both for the primary land-use models (H1 fixed effect)
# and for the substrate hypothesis (H4: CPOM vs k_total).

kinzig_micro <- kinzig %>% filter(type == "fine"   & is.finite(k) & k > 0)
kinzig_total <- kinzig %>% filter(type == "coarse" & is.finite(k) & k > 0)

# Pair coarse and fine bags by stream, replicate and land use to compute
# the macroinvertebrate contribution (k_macro = k_total - k_micro).
# Cases where the difference is non-positive are removed, as negative
# contributions are biologically implausible.

kinzig_macro <- kinzig_micro %>%
  select(stream, replicate, LandUse, k_microbial = k) %>%
  left_join(kinzig_total %>% select(stream, replicate, LandUse, k_total = k),
            by = c("stream", "replicate", "LandUse")) %>%
  mutate(k_macro = k_total - k_microbial) %>%
  filter(is.finite(k_macro), k_macro > 0)

message(glue("Kinzig – fine: {nrow(kinzig_micro)}, coarse: {nrow(kinzig_total)}, paired macro: {nrow(kinzig_macro)}"))


# 5 Boye catchment -------------------------------------------------------------

## 5.1 Restoration status assignment -------------------------------------------

# Source: SFB_BasicData_SamplingSites_Boye_&_Kinzig file.
# BOYkl     = Kleine Boye (reference, never restored)
# BOYohB224, BOYohKi, BOYuhHa = restored 2021 (~1 yr recovery at time of study)
# BOYohBr,   BOYohSp, BOYuhSp = restored 2002 (~20 yr recovery at time of study)

restoration_status <- tibble(
  stream = c("BOYkl", "BOYohB224", "BOYohBr", "BOYohKi",
             "BOYohSp", "BOYuhHa", "BOYuhSp"),
  RestStatus = c("reference",
                 "recently_restored",
                 "long_restored",
                 "recently_restored",
                 "long_restored",
                 "recently_restored",
                 "long_restored")
) %>%
  mutate(RestStatus = factor(RestStatus,
                             levels = c("reference", "recently_restored", "long_restored")))


## 5.2 Degree days and k – Boye ------------------------------------------------

# The Boye physico-chemistry file uses a different column structure from Kinzig.
# We rename columns to a common format before computing degree days.
# Verify column names with names(phch_boye) if the join fails.

phch_boye_summ <- phch_boye %>%
  rename(stream     = site,
         water_temp = water_temp_deg_c) %>%
  mutate(stream = str_trim(stream)) %>%
  group_by(stream) %>%
  summarise(water_temp = mean(water_temp, na.rm = TRUE), .groups = "drop")

boye <- processed %>%
  filter(catchment == "Boye") %>%
  left_join(restoration_status, by = "stream") %>%
  left_join(phch_boye_summ,     by = "stream") %>%
  mutate(
    deg_days = deployment_duration_days * water_temp,
    k        = -log(AFDM_end / AFDM_start) / deg_days
  ) %>%
  filter(is.finite(k), k > 0, AFDM_end < AFDM_start, deg_days > 0,
         !is.na(RestStatus)) %>%
  mutate(k = round(k, 5))

message(glue("Boye: {nrow(boye)} bags after filtering."))


## 5.3 Join master environmental predictors – Boye ------------------------------
#
# Bring in agricultural land use, conductivity, metal concentrations (for TU
# calculation in H2), and fine-substrate fractions (H3) from master_boye.
# We also include FPOM and sand (psam_pct) so we can classify sites by
# fine-substrate cover (>= 25% threshold, H3).

boye <- boye %>%
  left_join(
    master_boye %>%
      select(stream,
             pct_agri,
             conductivity_u_s_cm,
             no3_mg_l,
             ortho_po4_mg_l,
             total_po4_mg_l,
             cd_mg_l, cu_mg_l, fe_mg_l,
             ni_mg_l, pb_mg_l, zn_mg_l,
             fpom_pct,
             psam_pct),
    by = "stream"
  )


## 5.4 Prepare Boye analysis datasets -------------------------------------------
#
# Separate microbial (fine mesh) and total (coarse mesh) datasets.
# These feed the primary restoration-status models (Sections 7) and
# the substrate hypothesis (H3: fine sediment vs k_total).

boye_micro <- boye %>% filter(type == "fine"   & is.finite(k) & k > 0)
boye_total <- boye %>% filter(type == "coarse" & is.finite(k) & k > 0)

# Pair fine and coarse bags for macroinvertebrate contribution
boye_macro <- boye_micro %>%
  select(stream, replicate, RestStatus, k_microbial = k) %>%
  left_join(boye_total %>% select(stream, replicate, RestStatus, k_total = k),
            by = c("stream", "replicate", "RestStatus")) %>%
  mutate(k_macro = k_total - k_microbial) %>%
  filter(is.finite(k_macro), k_macro > 0)

message(glue("Boye – fine: {nrow(boye_micro)}, coarse: {nrow(boye_total)}, paired macro: {nrow(boye_macro)}"))


# =============================================================================
# 5a DLM calculation (Schreiner-style normalised mass loss)
# =============================================================================
message("=== Section 28: DLM calculation ===")

# DLM = % AFDM loss / (avg_temp_deg_c × deployment_duration_days)
# avg_temp_deg_c from master_vars (iButton logger).
# Temperature is encoded in the response → NOT used as a predictor.

temp_logger <- master_vars %>%
  select(stream, avg_temp_deg_c) %>%
  filter(is.finite(avg_temp_deg_c))

processed_dlm <- processed %>%
  left_join(temp_logger, by = "stream") %>%
  filter(
    is.finite(AFDM_start), AFDM_start > 0,
    is.finite(AFDM_end),   AFDM_end < AFDM_start,
    is.finite(avg_temp_deg_c), avg_temp_deg_c > 0,
    is.finite(deployment_duration_days), deployment_duration_days > 0
  ) %>%
  mutate(
    pct_mass_loss = (AFDM_start - AFDM_end) / AFDM_start * 100,
    DLM           = pct_mass_loss / (avg_temp_deg_c * deployment_duration_days)
  ) %>%
  filter(is.finite(DLM), DLM > 0)

message(glue("processed_dlm: {nrow(processed_dlm)} rows with valid DLM"))

## 28.1 Kinzig DLM datasets ---------------------------------------------------
dlm_kinzig_pred_cols <- c("pct_agri", "conductivity_u_s_cm", "no3_mg_l",
                           "ortho_po4_mg_l", "cpom_pct", "fpom_pct", "psam_pct",
                           "cd_mg_l", "cu_mg_l", "fe_mg_l",
                           "ni_mg_l", "pb_mg_l", "zn_mg_l")

kinzig_micro_dlm <- processed_dlm %>%
  filter(catchment == "Kinzig", type == "fine") %>%
  left_join(master_kinzig %>%
              select(stream, any_of(dlm_kinzig_pred_cols)) %>% distinct(),
            by = "stream") %>%
  left_join(site_info_kinzig, by = "stream") %>%
  filter(!is.na(LandUse)) %>%
  mutate(LandUse    = factor(LandUse, levels = c("rural", "urban")),
         fine_cover = replace_na(fpom_pct, 0) + replace_na(psam_pct, 0),
         fine_class = factor(if_else(fine_cover >= 25, "high_fine", "low_fine"),
                              levels = c("low_fine", "high_fine")))

kinzig_total_dlm <- processed_dlm %>%
  filter(catchment == "Kinzig", type == "coarse") %>%
  left_join(master_kinzig %>%
              select(stream, any_of(dlm_kinzig_pred_cols)) %>% distinct(),
            by = "stream") %>%
  left_join(site_info_kinzig, by = "stream") %>%
  filter(!is.na(LandUse)) %>%
  mutate(LandUse    = factor(LandUse, levels = c("rural", "urban")),
         fine_cover = replace_na(fpom_pct, 0) + replace_na(psam_pct, 0),
         fine_class = factor(if_else(fine_cover >= 25, "high_fine", "low_fine"),
                              levels = c("low_fine", "high_fine")))

kinzig_macro_dlm <- kinzig_micro_dlm %>%
  select(stream, replicate, LandUse, DLM_micro = DLM,
         any_of(dlm_kinzig_pred_cols)) %>%
  left_join(kinzig_total_dlm %>%
              select(stream, replicate, LandUse, DLM_total = DLM),
            by = c("stream", "replicate", "LandUse")) %>%
  mutate(DLM = DLM_total - DLM_micro) %>%
  filter(is.finite(DLM), DLM > 0)

## 28.2 Boye DLM datasets -----------------------------------------------------
dlm_boye_pred_cols <- c("pct_agri", "conductivity_u_s_cm", "no3_mg_l",
                         "ortho_po4_mg_l", "fpom_pct", "psam_pct",
                         "cd_mg_l", "cu_mg_l", "fe_mg_l",
                         "ni_mg_l", "pb_mg_l", "zn_mg_l")

pnec_dlm <- c(cd = 0.00025, cu = 0.0028, fe = 0.1,
              ni = 0.004,   pb = 0.0014, zn = 0.0078)

boye_micro_dlm <- processed_dlm %>%
  filter(catchment == "Boye", type == "fine") %>%
  left_join(master_boye %>%
              select(stream, any_of(dlm_boye_pred_cols)) %>% distinct(),
            by = "stream") %>%
  left_join(restoration_status, by = "stream") %>%
  filter(!is.na(RestStatus)) %>%
  mutate(
    fine_cover = replace_na(fpom_pct, 0) + replace_na(psam_pct, 0),
    fine_class = factor(if_else(fine_cover >= 25, "high_fine", "low_fine"),
                         levels = c("low_fine", "high_fine")),
    sumTU = rowSums(cbind(
      cd_mg_l / pnec_dlm["cd"], cu_mg_l / pnec_dlm["cu"],
      fe_mg_l / pnec_dlm["fe"], ni_mg_l / pnec_dlm["ni"],
      pb_mg_l / pnec_dlm["pb"], zn_mg_l / pnec_dlm["zn"]),
      na.rm = TRUE))

boye_total_dlm <- processed_dlm %>%
  filter(catchment == "Boye", type == "coarse") %>%
  left_join(master_boye %>%
              select(stream, any_of(dlm_boye_pred_cols)) %>% distinct(),
            by = "stream") %>%
  left_join(restoration_status, by = "stream") %>%
  filter(!is.na(RestStatus)) %>%
  mutate(
    fine_cover = replace_na(fpom_pct, 0) + replace_na(psam_pct, 0),
    fine_class = factor(if_else(fine_cover >= 25, "high_fine", "low_fine"),
                         levels = c("low_fine", "high_fine")),
    sumTU = rowSums(cbind(
      cd_mg_l / pnec_dlm["cd"], cu_mg_l / pnec_dlm["cu"],
      fe_mg_l / pnec_dlm["fe"], ni_mg_l / pnec_dlm["ni"],
      pb_mg_l / pnec_dlm["pb"], zn_mg_l / pnec_dlm["zn"]),
      na.rm = TRUE))

message(glue(
  "DLM datasets — Kinzig micro: {nrow(kinzig_micro_dlm)}, ",
  "total: {nrow(kinzig_total_dlm)}, macro: {nrow(kinzig_macro_dlm)} | ",
  "Boye micro: {nrow(boye_micro_dlm)}, total: {nrow(boye_total_dlm)}"
))




# 6 Mixed-effects models – Kinzig (land use effect) ----------------------------

# We model k as a function of LandUse with a random intercept for stream.
# The coefficient for LandUseurban gives the urban - rural difference.
# A negative estimate implies decomposition is lower in urban streams.

## 6.1 Fit models --------------------------------------------------------------

model_kinzig_macro <- lmer(k_macro     ~ LandUse + (1 | stream), data = kinzig_macro)
model_kinzig_micro <- lmer(k           ~ LandUse + (1 | stream), data = kinzig_micro)
model_kinzig_total <- lmer(k           ~ LandUse + (1 | stream), data = kinzig_total)


## 6.2 Model diagnostics -------------------------------------------------------

# Inspect residuals, QQ-plots etc. interactively before interpreting results.

check_kinzig_macro <- performance::check_model(model_kinzig_macro)
check_kinzig_micro <- performance::check_model(model_kinzig_micro)
check_kinzig_total <- performance::check_model(model_kinzig_total)


## 6.3 Extract and report fixed effects ----------------------------------------

tidy_kinzig_macro <- tidy(model_kinzig_macro, effects = "fixed")
tidy_kinzig_micro <- tidy(model_kinzig_micro, effects = "fixed")
tidy_kinzig_total <- tidy(model_kinzig_total, effects = "fixed")

extract_effect <- function(tdf, term_name = "LandUseurban") {
  tibble(
    est = tdf %>% filter(term == term_name) %>% pull(estimate),
    se  = tdf %>% filter(term == term_name) %>% pull(std.error),
    p   = tdf %>% filter(term == term_name) %>% pull(p.value)
  )
}

eff_kn_macro <- extract_effect(tidy_kinzig_macro)
eff_kn_micro <- extract_effect(tidy_kinzig_micro)
eff_kn_total <- extract_effect(tidy_kinzig_total)

message(glue("Kinzig H1 (k_macro): urban - rural = {round(eff_kn_macro$est, 5)} +/- {round(eff_kn_macro$se, 5)}; p = {signif(eff_kn_macro$p, 3)}"))
message(glue("Kinzig H2 (k_micro): urban - rural = {round(eff_kn_micro$est, 5)} +/- {round(eff_kn_micro$se, 5)}; p = {signif(eff_kn_micro$p, 3)}"))
message(glue("Kinzig H3 (k_total): urban - rural = {round(eff_kn_total$est, 5)} +/- {round(eff_kn_total$se, 5)}; p = {signif(eff_kn_total$p, 3)}"))


# 7 Mixed-effects models – Boye (restoration status effect) --------------------

# We model k as a function of RestStatus with a random intercept for stream.
# "reference" is the baseline level. Because RestStatus has three levels, we
# additionally compute pairwise contrasts with emmeans (Tukey-adjusted).

## 7.1 Fit models --------------------------------------------------------------

model_boye_macro <- lmer(k_macro ~ RestStatus + (1 | stream), data = boye_macro)
model_boye_micro <- lmer(k       ~ RestStatus + (1 | stream), data = boye_micro)
model_boye_total <- lmer(k       ~ RestStatus + (1 | stream), data = boye_total)


## 7.2 Model diagnostics -------------------------------------------------------

check_boye_macro <- performance::check_model(model_boye_macro)
check_boye_micro <- performance::check_model(model_boye_micro)
check_boye_total <- performance::check_model(model_boye_total)


## 7.3 Extract fixed effects ---------------------------------------------------

tidy_boye_macro <- tidy(model_boye_macro, effects = "fixed")
tidy_boye_micro <- tidy(model_boye_micro, effects = "fixed")
tidy_boye_total <- tidy(model_boye_total, effects = "fixed")


## 7.4 Pairwise contrasts (Tukey-adjusted) -------------------------------------

# Because RestStatus has three levels, we use emmeans to test all pairwise
# differences between reference, recently restored and long restored.

cont_boye_macro <- emmeans(model_boye_macro, pairwise ~ RestStatus)$contrasts %>%
  as.data.frame() %>% mutate(Endpoint = "k_macro")
cont_boye_micro <- emmeans(model_boye_micro, pairwise ~ RestStatus)$contrasts %>%
  as.data.frame() %>% mutate(Endpoint = "k_micro")
cont_boye_total <- emmeans(model_boye_total, pairwise ~ RestStatus)$contrasts %>%
  as.data.frame() %>% mutate(Endpoint = "k_total")

table_boye_contrasts <- bind_rows(cont_boye_macro, cont_boye_micro, cont_boye_total) %>%
  mutate(across(c(estimate, SE, t.ratio, p.value), ~ round(.x, 4)))

print(table_boye_contrasts)


# 8 Figures --------------------------------------------------------------------

## 8.1 Figure 1: Kinzig decomposition rates by land use ------------------------

p_kn_micro <- ggplot(kinzig_micro, aes(x = LandUse, y = k, fill = LandUse)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  ggtitle("") +
  ylab("Microbial k (per degree day)") +
  xlab("Land use") +
  labs(subtitle = "A") +
  theme_classic() +
  theme(axis.text.x    = element_text(colour = "black", size = 12),
        axis.text.y    = element_text(colour = "black", size = 12),
        plot.subtitle  = element_text(colour = "black", size = 14),
        legend.position = "none",
        text           = element_text(size = 11))

p_kn_total <- ggplot(kinzig_total, aes(x = LandUse, y = k, fill = LandUse)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  ggtitle("") +
  ylab("Total k (per degree day)") +
  xlab("Land use") +
  labs(subtitle = "B") +
  theme_classic() +
  theme(axis.text.x    = element_text(colour = "black", size = 12),
        axis.text.y    = element_text(colour = "black", size = 12),
        plot.subtitle  = element_text(colour = "black", size = 14),
        legend.position = "none",
        text           = element_text(size = 11))

p_kn_macro <- ggplot(kinzig_macro, aes(x = LandUse, y = k_macro, fill = LandUse)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  ggtitle("") +
  ylab(expression(k[macro]~"(per degree day)")) +
  xlab("Land use") +
  labs(subtitle = "C") +
  theme_classic() +
  theme(axis.text.x    = element_text(colour = "black", size = 12),
        axis.text.y    = element_text(colour = "black", size = 12),
        plot.subtitle  = element_text(colour = "black", size = 14),
        legend.position = "none",
        text           = element_text(size = 11))

png(filename = "Figure_1_Kinzig_decomposition_rates.png",
    width = 9, height = 3.5, units = "in", res = 900)
ggdraw() +
  draw_plot(p_kn_micro, x = 0,    y = 0, width = 0.33, height = 1) +
  draw_plot(p_kn_total, x = 0.33, y = 0, width = 0.33, height = 1) +
  draw_plot(p_kn_macro, x = 0.66, y = 0, width = 0.33, height = 1)
dev.off()


## 8.2 Figure 2: Boye decomposition rates by restoration status ----------------

rest_colours <- c("reference"         = "gray50",
                  "recently_restored" = "steelblue",
                  "long_restored"     = "forestgreen")

p_by_micro <- ggplot(boye_micro, aes(x = RestStatus, y = k, fill = RestStatus)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = rest_colours) +
  ggtitle("") +
  ylab("Microbial k (per degree day)") +
  xlab("Restoration status") +
  labs(subtitle = "A") +
  theme_classic() +
  theme(axis.text.x    = element_text(colour = "black", size = 10, angle = 15, hjust = 1),
        axis.text.y    = element_text(colour = "black", size = 12),
        plot.subtitle  = element_text(colour = "black", size = 14),
        legend.position = "none",
        text           = element_text(size = 11))

p_by_total <- ggplot(boye_total, aes(x = RestStatus, y = k, fill = RestStatus)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = rest_colours) +
  ggtitle("") +
  ylab("Total k (per degree day)") +
  xlab("Restoration status") +
  labs(subtitle = "B") +
  theme_classic() +
  theme(axis.text.x    = element_text(colour = "black", size = 10, angle = 15, hjust = 1),
        axis.text.y    = element_text(colour = "black", size = 12),
        plot.subtitle  = element_text(colour = "black", size = 14),
        legend.position = "none",
        text           = element_text(size = 11))

p_by_macro <- ggplot(boye_macro, aes(x = RestStatus, y = k_macro, fill = RestStatus)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = rest_colours) +
  ggtitle("") +
  ylab(expression(k[macro]~"(per degree day)")) +
  xlab("Restoration status") +
  labs(subtitle = "C") +
  theme_classic() +
  theme(axis.text.x    = element_text(colour = "black", size = 10, angle = 15, hjust = 1),
        axis.text.y    = element_text(colour = "black", size = 12),
        plot.subtitle  = element_text(colour = "black", size = 14),
        legend.position = "none",
        text           = element_text(size = 11))

png(filename = "Figure_2_Boye_decomposition_rates.png",
    width = 9, height = 3.5, units = "in", res = 900)
ggdraw() +
  draw_plot(p_by_micro, x = 0,    y = 0, width = 0.33, height = 1) +
  draw_plot(p_by_total, x = 0.33, y = 0, width = 0.33, height = 1) +
  draw_plot(p_by_macro, x = 0.66, y = 0, width = 0.33, height = 1)
dev.off()


# 9 Summary tables -------------------------------------------------------------

## 9.1 Mean decomposition rates by land use – Kinzig ---------------------------

summary_kinzig <- bind_rows(
  kinzig_micro %>%
    group_by(LandUse) %>%
    summarise(Catchment = "Kinzig", Endpoint = "Microbial", n = n(),
              mean_k = mean(k, na.rm = TRUE),
              se_k   = sd(k, na.rm = TRUE) / sqrt(n()), .groups = "drop"),
  kinzig_total %>%
    group_by(LandUse) %>%
    summarise(Catchment = "Kinzig", Endpoint = "Total", n = n(),
              mean_k = mean(k, na.rm = TRUE),
              se_k   = sd(k, na.rm = TRUE) / sqrt(n()), .groups = "drop"),
  kinzig_macro %>%
    group_by(LandUse) %>%
    summarise(Catchment = "Kinzig", Endpoint = "Macroinvertebrate", n = n(),
              mean_k = mean(k_macro, na.rm = TRUE),
              se_k   = sd(k_macro, na.rm = TRUE) / sqrt(n()), .groups = "drop")
)

print(summary_kinzig)


## 9.2 Mean decomposition rates by restoration status – Boye -------------------

summary_boye <- bind_rows(
  boye_micro %>%
    group_by(RestStatus) %>%
    summarise(Catchment = "Boye", Endpoint = "Microbial", n = n(),
              mean_k = mean(k, na.rm = TRUE),
              se_k   = sd(k, na.rm = TRUE) / sqrt(n()), .groups = "drop"),
  boye_total %>%
    group_by(RestStatus) %>%
    summarise(Catchment = "Boye", Endpoint = "Total", n = n(),
              mean_k = mean(k, na.rm = TRUE),
              se_k   = sd(k, na.rm = TRUE) / sqrt(n()), .groups = "drop"),
  boye_macro %>%
    group_by(RestStatus) %>%
    summarise(Catchment = "Boye", Endpoint = "Macroinvertebrate", n = n(),
              mean_k = mean(k_macro, na.rm = TRUE),
              se_k   = sd(k_macro, na.rm = TRUE) / sqrt(n()), .groups = "drop")
)

print(summary_boye)


## 9.3 Fixed effects tables ----------------------------------------------------

table_kinzig_models <- bind_rows(
  tidy_kinzig_macro %>% mutate(Catchment = "Kinzig", Endpoint = "Macroinvertebrate"),
  tidy_kinzig_micro %>% mutate(Catchment = "Kinzig", Endpoint = "Microbial"),
  tidy_kinzig_total %>% mutate(Catchment = "Kinzig", Endpoint = "Total")
)

table_boye_models <- bind_rows(
  tidy_boye_macro %>% mutate(Catchment = "Boye", Endpoint = "Macroinvertebrate"),
  tidy_boye_micro %>% mutate(Catchment = "Boye", Endpoint = "Microbial"),
  tidy_boye_total %>% mutate(Catchment = "Boye", Endpoint = "Total")
)

print(table_kinzig_models)
print(table_boye_models)


## 9.4 Save outputs to file ----------------------------------------------------

write.csv(summary_kinzig,         "Kinzig_decomposition_summary.csv",       row.names = FALSE)
write.csv(summary_boye,           "Boye_decomposition_summary.csv",         row.names = FALSE)
write.csv(table_kinzig_models,    "Kinzig_mixed_model_fixed_effects.csv",   row.names = FALSE)
write.csv(table_boye_models,      "Boye_mixed_model_fixed_effects.csv",     row.names = FALSE)
write.csv(table_boye_contrasts,   "Boye_pairwise_contrasts.csv",            row.names = FALSE)

ggsave("Kinzig_microbial_boxplot.png",  p_kn_micro, width = 3, height = 4, dpi = 300)
ggsave("Kinzig_total_boxplot.png",      p_kn_total, width = 3, height = 4, dpi = 300)
ggsave("Kinzig_macro_boxplot.png",      p_kn_macro, width = 3, height = 4, dpi = 300)
ggsave("Boye_microbial_boxplot.png",    p_by_micro, width = 3, height = 4, dpi = 300)
ggsave("Boye_total_boxplot.png",        p_by_total, width = 3, height = 4, dpi = 300)
ggsave("Boye_macro_boxplot.png",        p_by_macro, width = 3, height = 4, dpi = 300)

# =============================================================================
# 10  H1 – Agricultural land use and conductivity across both catchments
# =============================================================================
#
# We test whether the percentage of agricultural land in the upstream catchment
# and conductivity (an integrative measure of ion loading) jointly predict
# microbial decomposition rate across all study streams. Both predictors are
# scaled to zero mean and unit variance (scale()) so that their coefficients
# are directly comparable. The random intercept for stream accounts for the
# fact that multiple leaf bags were deployed per site.
#
# A negative estimate for pct_agri or conductivity would support H1.
# Reference: Schreiner et al. (2023); Bundschuh et al. (2021).
# =============================================================================

## 10.1 Build cross-catchment microbial dataset --------------------------------
#
# We pool boye_micro and kinzig_micro, keeping only the columns needed for the
# model. Both datasets now carry pct_agri and conductivity_u_s_cm after the
# joins in sections 4.3 and 5.3.

all_micro_H1 <- bind_rows(
  boye_micro   %>% select(stream, catchment, k, pct_agri, conductivity_u_s_cm),
  kinzig_micro %>% select(stream, catchment, k, pct_agri, conductivity_u_s_cm)
) %>%
  filter(is.finite(k), k > 0,
         is.finite(pct_agri),
         is.finite(conductivity_u_s_cm)) %>%
  mutate(
    pct_agri_sc    = scale(pct_agri)[, 1],
    conductivity_sc = scale(conductivity_u_s_cm)[, 1]
  )

message(glue("H1 dataset: {nrow(all_micro_H1)} bags across {n_distinct(all_micro_H1$stream)} streams"))


## 10.2 Fit H1 model -----------------------------------------------------------

model_H1 <- lmer(k ~ pct_agri_sc + conductivity_sc + (1 | stream),
                 data = all_micro_H1)

summary(model_H1)
check_H1 <- performance::check_model(model_H1)
tidy_H1  <- tidy(model_H1, effects = "fixed")
print(tidy_H1)


## 10.3 H1 figure – partial effects plot ---------------------------------------
#
# ggeffects::ggpredict() computes the marginal effect of each predictor while
# holding the other at its mean. We plot both effects side by side.

eff_H1_agri <- ggpredict(model_H1, terms = "pct_agri_sc [all]")
eff_H1_cond <- ggpredict(model_H1, terms = "conductivity_sc [all]")

p_H1_agri <- plot(eff_H1_agri) +
  xlab("Agricultural land use (scaled)") +
  ylab("Predicted microbial k") +
  ggtitle("H1: Agriculture effect") +
  theme_classic()

p_H1_cond <- plot(eff_H1_cond) +
  xlab("Conductivity (scaled)") +
  ylab("Predicted microbial k") +
  ggtitle("H1: Conductivity effect") +
  theme_classic()

png("Figure_H1_agriculture_conductivity.png", width = 8, height = 4,
    units = "in", res = 300)
plot_grid(p_H1_agri, p_H1_cond, labels = c("A", "B"), nrow = 1)
dev.off()

write.csv(tidy_H1, "H1_model_fixed_effects.csv", row.names = FALSE)


# =============================================================================
# 11  H2 – Metal toxic units (Boye) vs nutrient concentrations (Kinzig)
# =============================================================================
#
# We test the catchment-specific stressor hypotheses:
#   Boye:   metal contamination (sum of toxic units) reduces k_microbial.
#   Kinzig: elevated nutrients (NO3, ortho-PO4) enhance k_microbial.
#
# Toxic units are calculated as:
#   TU_metal = concentration_measured / PNEC_metal
# where PNEC values are EU Environmental Quality Standards (annual average,
# freshwater, dissolved phase):
#   Cd = 0.00025 mg/L, Cu = 0.0028 mg/L, Fe = 0.1 mg/L (Dutch guideline),
#   Ni = 0.004  mg/L,  Pb = 0.0014 mg/L, Zn = 0.0078 mg/L
#
# A negative TU estimate supports H2 for Boye; positive NO3/PO4 estimates
# support H2 for Kinzig.
# =============================================================================

## 11.1 Boye – compute metal TU and fit model ----------------------------------
#
# Each metal TU is concentration / PNEC. The sumTU aggregates overall toxicity
# pressure, which is the standard approach when multiple metals co-occur.

pnec <- c(cd = 0.00025, cu = 0.0028, fe = 0.1,
          ni = 0.004,   pb = 0.0014, zn = 0.0078)

boye_micro_H2 <- boye_micro %>%
  filter(is.finite(cd_mg_l), is.finite(cu_mg_l), is.finite(fe_mg_l),
         is.finite(ni_mg_l), is.finite(pb_mg_l), is.finite(zn_mg_l)) %>%
  mutate(
    TU_cd  = cd_mg_l / pnec["cd"],
    TU_cu  = cu_mg_l / pnec["cu"],
    TU_fe  = fe_mg_l / pnec["fe"],
    TU_ni  = ni_mg_l / pnec["ni"],
    TU_pb  = pb_mg_l / pnec["pb"],
    TU_zn  = zn_mg_l / pnec["zn"],
    sumTU  = TU_cd + TU_cu + TU_fe + TU_ni + TU_pb + TU_zn,
    sumTU_sc = scale(sumTU)[, 1]
  ) %>%
  filter(is.finite(sumTU), sumTU > 0)

message(glue("H2 Boye metal TU dataset: {nrow(boye_micro_H2)} bags, {n_distinct(boye_micro_H2$stream)} streams"))

model_H2_boye <- lmer(k ~ sumTU_sc + (1 | stream), data = boye_micro_H2)
summary(model_H2_boye)
tidy_H2_boye <- tidy(model_H2_boye, effects = "fixed")
print(tidy_H2_boye)


## 11.2 Kinzig – nutrient model ------------------------------------------------
#
# We test NO3 and ortho-PO4 as separate predictors. Both are scaled.
# A positive estimate for either nutrient supports the enrichment-stimulation
# hypothesis for fungal and bacterial decomposers.

kinzig_micro_H2 <- kinzig_micro %>%
  filter(is.finite(no3_mg_l), is.finite(ortho_po4_mg_l)) %>%
  mutate(
    no3_sc  = scale(no3_mg_l)[, 1],
    po4_sc  = scale(ortho_po4_mg_l)[, 1]
  )

message(glue("H2 Kinzig nutrient dataset: {nrow(kinzig_micro_H2)} bags, {n_distinct(kinzig_micro_H2$stream)} streams"))

model_H2_kinzig <- lmer(k ~ no3_sc + po4_sc + (1 | stream),
                        data = kinzig_micro_H2)
summary(model_H2_kinzig)
tidy_H2_kinzig <- tidy(model_H2_kinzig, effects = "fixed")
print(tidy_H2_kinzig)


## 11.3 H2 figures -------------------------------------------------------------

eff_H2_boye   <- ggpredict(model_H2_boye,   terms = "sumTU_sc [all]")
eff_H2_no3    <- ggpredict(model_H2_kinzig,  terms = "no3_sc [all]")
eff_H2_po4    <- ggpredict(model_H2_kinzig,  terms = "po4_sc [all]")

p_H2_boye <- plot(eff_H2_boye) +
  xlab("Sum metal TU (scaled)") +
  ylab("Predicted microbial k") +
  ggtitle("H2 Boye: metal toxicity") +
  theme_classic()

p_H2_no3 <- plot(eff_H2_no3) +
  xlab(expression(NO[3]^{"-"}~"(scaled)")) +
  ylab("Predicted microbial k") +
  ggtitle(expression("H2 Kinzig: NO"[3])) +
  theme_classic()

p_H2_po4 <- plot(eff_H2_po4) +
  xlab("Ortho-PO4 (scaled)") +
  ylab("Predicted microbial k") +
  ggtitle("H2 Kinzig: Ortho-PO4") +
  theme_classic()

png("Figure_H2_metals_nutrients.png", width = 10, height = 4,
    units = "in", res = 300)
plot_grid(p_H2_boye, p_H2_no3, p_H2_po4, labels = c("A", "B", "C"), nrow = 1)
dev.off()

write.csv(tidy_H2_boye,   "H2_Boye_metal_TU_model.csv",    row.names = FALSE)
write.csv(tidy_H2_kinzig, "H2_Kinzig_nutrient_model.csv",  row.names = FALSE)


# =============================================================================
# 12  H3 – Fine substrate and total decomposition in the Boye
# =============================================================================
#
# We classify each Boye site as "high_fine" (FPOM + sand >= 25% streambed
# cover) or "low_fine" (< 25%). The 25% threshold follows Lepori et al.
# (2005), who showed that macroinvertebrate shredder communities are
# substantially impaired above this level of fine sediment cover.
#
# We compare k_total between the two groups using a linear mixed-effects model
# with a random intercept for stream. A negative estimate for high_fine would
# support H3.
# =============================================================================

## 12.1 Classify fine substrate -------------------------------------------------

boye_total_H3 <- boye_total %>%
  mutate(
    fine_cover  = replace_na(fpom_pct, 0) + replace_na(psam_pct, 0),
    fine_class  = factor(
      ifelse(fine_cover >= 25, "high_fine", "low_fine"),
      levels = c("low_fine", "high_fine")   # low_fine = reference level
    )
  ) %>%
  filter(is.finite(fine_cover))

message(glue(
  "H3 fine-substrate groups: ",
  "low_fine = {sum(boye_total_H3$fine_class == 'low_fine')}, ",
  "high_fine = {sum(boye_total_H3$fine_class == 'high_fine')}"
))


## 12.2 Fit H3 model ------------------------------------------------------------

model_H3 <- lmer(k ~ fine_class + (1 | stream), data = boye_total_H3)
summary(model_H3)
tidy_H3 <- tidy(model_H3, effects = "fixed")
print(tidy_H3)


## 12.3 H3 figure ---------------------------------------------------------------

p_H3 <- ggplot(boye_total_H3, aes(x = fine_class, y = k, fill = fine_class)) +
  geom_boxplot(width = 0.5, alpha = 0.85) +
  scale_fill_manual(values = c("low_fine"  = "steelblue",
                               "high_fine" = "sandybrown")) +
  xlab("Fine substrate class") +
  ylab("Total k (per degree day)") +
  ggtitle("H3: Fine substrate vs total decomposition (Boye)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        text      = element_text(size = 11))

png("Figure_H3_fine_substrate_Boye.png", width = 5, height = 4,
    units = "in", res = 300)
print(p_H3)
dev.off()

write.csv(tidy_H3, "H3_fine_substrate_model.csv", row.names = FALSE)


# =============================================================================
# 13  H4 – CPOM availability and total decomposition in the Kinzig
# =============================================================================
#
# We test whether the percentage of coarse particulate organic matter (CPOM)
# on the streambed positively predicts total leaf decomposition across Kinzig
# sites. CPOM provides both habitat and food for shredding macroinvertebrates,
# so a higher CPOM percentage is expected to support a more abundant and
# diverse shredder community and thus a higher k_total (Graca et al., 2001).
#
# We fit a continuous regression of k_total on cpom_pct, with stream as a
# random intercept to account for multiple replicates per site.
# =============================================================================

## 13.1 Prepare H4 dataset ------------------------------------------------------

kinzig_total_H4 <- kinzig_total %>%
  filter(is.finite(cpom_pct)) %>%
  mutate(cpom_sc = scale(cpom_pct)[, 1])

message(glue("H4 CPOM dataset: {nrow(kinzig_total_H4)} bags, {n_distinct(kinzig_total_H4$stream)} streams"))


## 13.2 Fit H4 model ------------------------------------------------------------

model_H4 <- lmer(k ~ cpom_sc + (1 | stream), data = kinzig_total_H4)
summary(model_H4)
tidy_H4 <- tidy(model_H4, effects = "fixed")
print(tidy_H4)


## 13.3 H4 figure ---------------------------------------------------------------
#
# Scatter plot of observed k_total against CPOM coverage, with the model
# partial-effect regression line drawn through the points.

eff_H4 <- ggpredict(model_H4, terms = "cpom_sc [all]")

p_H4 <- ggplot(kinzig_total_H4, aes(x = cpom_sc, y = k)) +
  geom_point(aes(colour = LandUse), size = 2.5, alpha = 0.8) +
  geom_ribbon(data = as.data.frame(eff_H4),
              aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, alpha = 0.2, fill = "gray40") +
  geom_line(data  = as.data.frame(eff_H4),
            aes(x = x, y = predicted),
            inherit.aes = FALSE, linewidth = 0.9, colour = "black") +
  scale_colour_manual(values = c("rural" = "forestgreen", "urban" = "darkorange")) +
  xlab("CPOM cover (scaled)") +
  ylab("Total k (per degree day)") +
  ggtitle("H4: CPOM availability vs total decomposition (Kinzig)") +
  theme_classic() +
  theme(legend.title = element_text(size = 10),
        axis.text    = element_text(size = 12),
        text         = element_text(size = 11))

png("Figure_H4_CPOM_Kinzig.png", width = 6, height = 4,
    units = "in", res = 300)
print(p_H4)
dev.off()

write.csv(tidy_H4, "H4_CPOM_model.csv", row.names = FALSE)


# =============================================================================
# 14  Schreiner Step 2 – Substrate comparison and response distribution check
# =============================================================================
#
# Following Schreiner et al., we first check whether the two substrates
# (fine mesh = microbial, coarse mesh = total) produce correlated responses.
# A high Pearson r (> 0.7) means they are interchangeable; a low r justifies
# analysing them separately. We also inspect the skewness of each k
# distribution: if |skewness| > 2 we apply a log10 transformation before
# variable selection and modelling, because LASSO and OLS assume approximate
# normality in the response. We additionally report inter-replicate variance
# per site as a data-quality diagnostic (high relative variance flags
# problematic sites).
# =============================================================================

## 14.1 Aggregate bag-level k to site-level means ------------------------------
#
# Predictors in master_vars are measured once per site, so we reduce the
# response to one value per site before merging. We do this separately for
# each catchment and endpoint.

site_k_kinzig <- kinzig_micro %>%
  group_by(stream) %>%
  summarise(k_micro = mean(k, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    kinzig_total %>%
      group_by(stream) %>%
      summarise(k_total = mean(k, na.rm = TRUE), .groups = "drop"),
    by = "stream"
  )

site_k_boye <- boye_micro %>%
  group_by(stream) %>%
  summarise(k_micro = mean(k, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    boye_total %>%
      group_by(stream) %>%
      summarise(k_total = mean(k, na.rm = TRUE), .groups = "drop"),
    by = "stream"
  )


## 14.2 Pearson correlation between substrates ---------------------------------

cor_kinzig <- cor.test(site_k_kinzig$k_micro, site_k_kinzig$k_total,
                       method = "pearson", use = "complete.obs")
cor_boye   <- cor.test(site_k_boye$k_micro,   site_k_boye$k_total,
                       method = "pearson", use = "complete.obs")

message(glue("Kinzig micro-total Pearson r = {round(cor_kinzig$estimate, 3)}, ",
             "p = {signif(cor_kinzig$p.value, 3)}"))
message(glue("Boye micro-total Pearson r = {round(cor_boye$estimate, 3)}, ",
             "p = {signif(cor_boye$p.value, 3)}"))


## 14.3 Skewness check and conditional log10 transformation --------------------
#
# We compute skewness at the bag level (full sample) for each response
# variable. A |skewness| > 2 triggers log10 transformation. The transformed
# value (k_resp) is stored back into the bag-level dataframes and used in all
# downstream modelling sections.

skew_table <- tibble(
  dataset  = c("kinzig_micro k", "kinzig_total k",
                "boye_micro k",   "boye_total k"),
  skewness = c(
    e1071::skewness(kinzig_micro$k, na.rm = TRUE),
    e1071::skewness(kinzig_total$k, na.rm = TRUE),
    e1071::skewness(boye_micro$k,   na.rm = TRUE),
    e1071::skewness(boye_total$k,   na.rm = TRUE)
  )
) %>%
  mutate(log_transform_applied = abs(skewness) > 2)

print(skew_table)

# Apply transformation and store as k_resp
add_k_resp <- function(df) {
  sk <- e1071::skewness(df$k, na.rm = TRUE)
  df %>% mutate(k_resp = if (abs(sk) > 2) log10(k) else k)
}

kinzig_micro <- add_k_resp(kinzig_micro)
kinzig_total <- add_k_resp(kinzig_total)
boye_micro   <- add_k_resp(boye_micro)
boye_total   <- add_k_resp(boye_total)


## 14.4 Inter-replicate variability per site -----------------------------------
#
# Relative variance (var / mean^2) close to zero indicates consistent
# replicates. Values > 0.5 flag sites where within-site variability is large
# relative to the mean, which can inflate residuals in models.

replicate_var_kinzig <- kinzig_micro %>%
  group_by(stream, LandUse) %>%
  summarise(
    n       = n(),
    mean_k  = mean(k, na.rm = TRUE),
    var_k   = var(k,  na.rm = TRUE),
    rel_var = var_k / mean_k^2,
    .groups = "drop"
  )

replicate_var_boye <- boye_micro %>%
  group_by(stream, RestStatus) %>%
  summarise(
    n       = n(),
    mean_k  = mean(k, na.rm = TRUE),
    var_k   = var(k,  na.rm = TRUE),
    rel_var = var_k / mean_k^2,
    .groups = "drop"
  )

print(replicate_var_kinzig)
print(replicate_var_boye)

write.csv(skew_table,            "Response_skewness_check.csv",       row.names = FALSE)
write.csv(replicate_var_kinzig,  "Kinzig_replicate_variability.csv",  row.names = FALSE)
write.csv(replicate_var_boye,    "Boye_replicate_variability.csv",    row.names = FALSE)


# =============================================================================
# 15  Schreiner Step 3 – Prepare explanatory variables
# =============================================================================
#
# We build one site-level predictor table per catchment by joining the
# site-mean k values from section 14 with the master_vars predictors. We then
# apply the following sequential filters:
#
#   (i)   Remove temperature variables (already embedded in k via degree days).
#   (ii)  Remove variables with near-zero range (max - min < 1e-6) or where
#         > 50% of values are exactly 0 – these carry no information for
#         regression or LASSO.
#   (iii) Check skewness of each numeric predictor; log10-transform those
#         with |skewness| > 2 (add 1 first if the variable contains zeros).
#   (iv)  Compute pairwise Pearson correlations; report all pairs with r > 0.7.
#         Where two variables are strongly correlated, retain the one with
#         greater ecological relevance to H1-H4.
#   (v)   Fit a full OLS model on k_micro, compute VIF; drop variables > 10.
# =============================================================================

## 15.1 Build site-level tables ------------------------------------------------
#
# Temperature columns are removed here because the decomposition rate k is
# already normalised by temperature x days; retaining temperature as a
# predictor would create a spurious negative relationship (higher T = lower k
# even when biology is unchanged).

excl_always <- c("stream", "catchment",
                 "water_temp_deg_c", "avg_temp_deg_c",
                 "comments", "date")

pred_cols_kinzig <- setdiff(names(master_kinzig), excl_always)
pred_cols_boye   <- setdiff(names(master_boye),   excl_always)

site_kinzig <- site_k_kinzig %>%
  left_join(master_kinzig %>% select(stream, all_of(pred_cols_kinzig)),
            by = "stream")

site_boye <- site_k_boye %>%
  left_join(master_boye %>% select(stream, all_of(pred_cols_boye)),
            by = "stream")

message(glue("Site-level table – Kinzig: {nrow(site_kinzig)} sites, ",
             "{ncol(site_kinzig) - 3} predictors before filtering"))
message(glue("Site-level table – Boye: {nrow(site_boye)} sites, ",
             "{ncol(site_boye) - 3} predictors before filtering"))


## 15.2 Remove low-information variables ---------------------------------------

filter_low_info <- function(df, resp_cols) {
  pred_df <- df %>% select(-stream, -all_of(resp_cols))
  keep    <- sapply(pred_df, function(x) {
    xn <- x[is.finite(x)]
    if (length(xn) == 0) return(FALSE)
    (max(xn) - min(xn)) > 1e-6 & mean(xn == 0, na.rm = TRUE) < 0.50
  })
  c("stream", resp_cols, names(pred_df)[keep])
}

resp_cols <- c("k_micro", "k_total")

site_kinzig <- site_kinzig %>% select(all_of(filter_low_info(site_kinzig, resp_cols)))
site_boye   <- site_boye   %>% select(all_of(filter_low_info(site_boye,   resp_cols)))

message(glue("After low-info filter – Kinzig: {ncol(site_kinzig) - 3} predictors remaining"))
message(glue("After low-info filter – Boye:   {ncol(site_boye)   - 3} predictors remaining"))


## 15.3 Skewness check and log10 transformation of predictors ------------------
#
# We transform each predictor with |skewness| > 2. Where the minimum value
# is 0 we use log10(x + 1) to avoid -Inf. This matches the procedure applied
# to the response variable in section 14.

transform_predictors <- function(df, skip_cols) {
  df_out <- df
  for (col in setdiff(names(df), skip_cols)) {
    x <- df[[col]]
    if (!is.numeric(x)) next
    sk <- e1071::skewness(x[is.finite(x)], na.rm = TRUE)
    if (!is.na(sk) && abs(sk) > 2) {
      if (min(x[is.finite(x)]) == 0) {
        df_out[[col]] <- log10(x + 1)
        message(glue("  log10(x+1) -> {col}  (skew = {round(sk, 2)})"))
      } else {
        df_out[[col]] <- log10(x)
        message(glue("  log10(x)   -> {col}  (skew = {round(sk, 2)})"))
      }
    }
  }
  df_out
}

skip_resp <- c("stream", "k_micro", "k_total")

message("Kinzig predictor transformations:")
site_kinzig <- transform_predictors(site_kinzig, skip_resp)

message("Boye predictor transformations:")
site_boye   <- transform_predictors(site_boye,   skip_resp)


## 15.4 Collinearity matrix (Pearson r > 0.7) ----------------------------------
#
# We flag all predictor pairs whose absolute Pearson correlation exceeds 0.7.
# For each flagged pair we retain the variable that is directly referenced
# in H1-H4 or that is more ecologically interpretable at catchment scale.
# The decisions are documented in the comments below and applied in 15.5.

flag_collinear <- function(df, skip_cols, threshold = 0.7) {
  num_df   <- df %>% select(-all_of(skip_cols)) %>%
    select(where(is.numeric)) %>% na.omit()
  corr_mat <- cor(num_df, method = "pearson")
  idx      <- which(abs(corr_mat) > threshold & lower.tri(corr_mat), arr.ind = TRUE)
  if (nrow(idx) == 0) { message("No collinear pairs detected."); return(tibble()) }
  tibble(
    var1 = rownames(corr_mat)[idx[, 1]],
    var2 = colnames(corr_mat)[idx[, 2]],
    r    = round(corr_mat[idx], 3)
  ) %>% arrange(desc(abs(r)))
}

message("Kinzig collinear predictor pairs (r > 0.7):")
collinear_kinzig <- flag_collinear(site_kinzig, skip_resp)
print(collinear_kinzig)

message("Boye collinear predictor pairs (r > 0.7):")
collinear_boye   <- flag_collinear(site_boye,   skip_resp)
print(collinear_boye)

write.csv(collinear_kinzig, "Kinzig_collinear_pairs.csv", row.names = FALSE)
write.csv(collinear_boye,   "Boye_collinear_pairs.csv",   row.names = FALSE)


## 15.5 Manual collinearity resolution ----------------------------------------
#
# Raw land-use areas (area_agric, area_urban, area_rural, area_indus) are
# redundant with their percentage equivalents (pct_agri, pct_urban, etc.),
# which are more comparable across catchments of different sizes -> drop raw.
#
# total_po4_mg_l correlates with ortho_po4_mg_l (reactive fraction is more
# directly bioavailable to decomposers) -> keep ortho_po4_mg_l.
#
# tn_mg_l correlates with no3_mg_l in agricultural streams (NO3 dominates
# the nitrogen load) -> keep no3_mg_l (direct H2 predictor).
#
# o2_sat_pct correlates with do_mg_l -> keep do_mg_l (absolute concentration).
#
# hco3_mg_l and co3_mg_l are carbonate system pairs; keep hco3_mg_l
# (bicarbonate, the dominant form in near-neutral freshwaters).

drop_collinear <- c("area_agric", "area_indus", "area_urban", "area_rural",
                    "total_po4_mg_l", "tn_mg_l", "o2_sat_pct", "co3_mg_l")

site_kinzig <- site_kinzig %>% select(-any_of(drop_collinear))
site_boye   <- site_boye   %>% select(-any_of(drop_collinear))

message(glue("After collinearity filter – Kinzig: {ncol(site_kinzig) - 3} predictors"))
message(glue("After collinearity filter – Boye:   {ncol(site_boye)   - 3} predictors"))


## 15.6 VIF check and removal of variables > 10 --------------------------------
#
# We fit a full OLS model on k_micro against all remaining numeric predictors,
# using complete cases only. VIF > 10 indicates severe multicollinearity that
# would destabilise LASSO coefficient paths. We drop flagged variables.
# Note: with only 7 Boye sites the full model may be rank-deficient; we print
# a warning if VIF cannot be computed.

compute_vif <- function(df, response, skip_cols) {
  df_cc <- df %>%
    select(-all_of(setdiff(skip_cols, response))) %>%
    select(where(is.numeric)) %>%
    na.omit()
  if (nrow(df_cc) <= ncol(df_cc)) {
    message(glue("  VIF skipped: n ({nrow(df_cc)}) <= p ({ncol(df_cc)}). ",
                 "Reduce predictors manually before proceeding."))
    return(tibble(variable = character(), VIF = numeric()))
  }
  fmla     <- as.formula(paste(response, "~ ."))
  lm_full  <- lm(fmla, data = df_cc)
  vif_vals <- car::vif(lm_full)
  tibble(variable = names(vif_vals), VIF = round(as.numeric(vif_vals), 2)) %>%
    arrange(desc(VIF))
}

message("Kinzig VIF (k_micro):")
vif_kinzig <- compute_vif(site_kinzig, "k_micro", c("stream", "k_total"))
print(vif_kinzig)

message("Boye VIF (k_micro):")
vif_boye <- compute_vif(site_boye, "k_micro", c("stream", "k_total"))
print(vif_boye)

# Drop VIF > 10
drop_vif_kinzig <- vif_kinzig %>% filter(VIF > 10) %>% pull(variable)
drop_vif_boye   <- vif_boye   %>% filter(VIF > 10) %>% pull(variable)

if (length(drop_vif_kinzig) > 0) {
  message(glue("  Dropping Kinzig VIF>10: {paste(drop_vif_kinzig, collapse = ', ')}"))
  site_kinzig <- site_kinzig %>% select(-any_of(drop_vif_kinzig))
}
if (length(drop_vif_boye) > 0) {
  message(glue("  Dropping Boye VIF>10: {paste(drop_vif_boye, collapse = ', ')}"))
  site_boye <- site_boye %>% select(-any_of(drop_vif_boye))
}

write.csv(vif_kinzig, "Kinzig_VIF_table.csv", row.names = FALSE)
write.csv(vif_boye,   "Boye_VIF_table.csv",   row.names = FALSE)


# =============================================================================
# 16  Schreiner Step 4 – Variable selection via LASSO stability selection
# =============================================================================
#
# Stability selection (Meinshausen & Buhlmann, 2010) runs LASSO repeatedly
# on bootstrapped subsamples and records how often each variable is selected
# at the optimal penalty. Variables selected in >= 70% of runs (cutoff = 0.70)
# with a per-family error rate (PFER) of 1 are considered robustly predictive.
#
# We run stability selection separately for k_micro and k_total in each
# catchment. All predictors are standardised to zero mean and unit variance
# before fitting so that LASSO penalises them on a common scale.
#
# Important note on sample size: the Boye dataset has only 7 sites. With
# n = 7, stability selection results should be treated as exploratory and
# interpreted alongside the biological reasoning in H1-H4.
# =============================================================================

run_stabsel <- function(df, response, skip_cols, cutoff = 0.70, PFER = 1) {
  df_cc <- df %>%
    select(-all_of(setdiff(skip_cols, response))) %>%
    select(where(is.numeric)) %>%
    na.omit()

  y <- df_cc[[response]]
  X <- df_cc %>% select(-all_of(response))

  # Remove zero-variance columns after na.omit
  X <- X %>% select(where(~ sd(.x, na.rm = TRUE) > 0))

  if (ncol(X) < 2 || nrow(X) < 5) {
    message(glue("  Stability selection skipped for {response}: ",
                 "n = {nrow(X)}, p = {ncol(X)} (insufficient data)."))
    return(NULL)
  }

  X_sc <- scale(X)

  stab <- stabsel(
    x          = X_sc,
    y          = y,
    fitfun     = glmnet.lasso,
    cutoff     = cutoff,
    PFER       = PFER,
    assumption = "unimodal"
  )
  stab
}

skip_base <- "stream"

message("Running stability selection – Kinzig k_micro ...")
stab_kn_micro <- run_stabsel(site_kinzig, "k_micro", c(skip_base, "k_total"))

message("Running stability selection – Kinzig k_total ...")
stab_kn_total <- run_stabsel(site_kinzig, "k_total", c(skip_base, "k_micro"))

message("Running stability selection – Boye k_micro ...")
stab_by_micro <- run_stabsel(site_boye, "k_micro", c(skip_base, "k_total"))

message("Running stability selection – Boye k_total ...")
stab_by_total <- run_stabsel(site_boye, "k_total", c(skip_base, "k_micro"))


## 16.1 Extract selected variables per response --------------------------------

extract_selected <- function(stab, label) {
  if (is.null(stab)) {
    return(tibble(dataset = label, variable = NA_character_,
                  selection_prob = NA_real_))
  }
  tibble(
    dataset        = label,
    variable       = names(stab$max),
    selection_prob = as.numeric(stab$max)
  ) %>%
    filter(selection_prob >= stab$cutoff) %>%
    arrange(desc(selection_prob))
}

selected_all <- bind_rows(
  extract_selected(stab_kn_micro, "Kinzig k_micro"),
  extract_selected(stab_kn_total, "Kinzig k_total"),
  extract_selected(stab_by_micro, "Boye k_micro"),
  extract_selected(stab_by_total, "Boye k_total")
)

message("Stability-selected variables:")
print(selected_all)
write.csv(selected_all, "Stability_selection_results.csv", row.names = FALSE)


## 16.2 Stability path plots ---------------------------------------------------
#
# The stability path shows selection probability vs LASSO penalty for each
# variable. Variables crossing the cutoff line are considered selected.

save_stabsel_plot <- function(stab, filename, main_title) {
  if (is.null(stab)) return(invisible(NULL))
  png(filename, width = 7, height = 5, units = "in", res = 300)
  plot(stab, main = main_title)
  dev.off()
  message(glue("Saved: {filename}"))
}

save_stabsel_plot(stab_kn_micro, "Stabsel_Kinzig_k_micro.png",
                  "Stability selection: Kinzig microbial k")
save_stabsel_plot(stab_kn_total, "Stabsel_Kinzig_k_total.png",
                  "Stability selection: Kinzig total k")
save_stabsel_plot(stab_by_micro, "Stabsel_Boye_k_micro.png",
                  "Stability selection: Boye microbial k")
save_stabsel_plot(stab_by_total, "Stabsel_Boye_k_total.png",
                  "Stability selection: Boye total k")


# =============================================================================
# 17  Schreiner Step 5 – Marginal effects of stability-selected variables
# =============================================================================
#
# For each catchment and response variable we fit a full OLS model including
# all stability-selected predictors, then plot the marginal effect of each
# selected variable using ggpredict() while holding all other predictors at
# their mean. We use lm() here (as in Schreiner et al.) because at the site
# level (one observation per stream) the mixed-model random intercept is no
# longer estimable. The resulting figures directly correspond to the
# hypothesis predictions in H1-H4 and are the primary output figures for the
# thesis Results section.
# =============================================================================

fit_and_plot_marginal <- function(df, response, selected_vars,
                                   skip_cols, catchment_label,
                                   response_label, response_axis) {
  # Guard: need at least one selected variable that exists in df
  avail <- intersect(selected_vars[!is.na(selected_vars)], names(df))
  if (length(avail) == 0) {
    message(glue("  No selected vars available for {catchment_label} {response_label}"))
    return(invisible(NULL))
  }

  df_cc  <- df %>%
    select(all_of(c(response, avail))) %>%
    na.omit()

  if (nrow(df_cc) < length(avail) + 1) {
    message(glue("  Not enough complete cases for {catchment_label} {response_label}"))
    return(invisible(NULL))
  }

  fmla   <- as.formula(paste(response, "~", paste(avail, collapse = " + ")))
  lm_fit <- lm(fmla, data = df_cc)

  message(glue("{catchment_label} {response_label} marginal effects model:"))
  print(summary(lm_fit))

  # One ggpredict panel per selected variable
  panels <- lapply(avail, function(v) {
    eff <- ggpredict(lm_fit, terms = glue("{v} [all]"))
    plot(eff) +
      xlab(gsub("_", " ", v)) +
      ylab(response_axis) +
      ggtitle(glue("{catchment_label}: {gsub('_', ' ', v)}")) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(size = 9))
  })

  ncols  <- min(3L, length(panels))
  height <- ceiling(length(panels) / ncols) * 3.5
  fig    <- plot_grid(plotlist = panels, ncol = ncols)

  fname  <- glue("Figure_Marginal_{catchment_label}_{response_label}.png")
  png(fname, width = ncols * 4, height = height, units = "in", res = 300)
  print(fig)
  dev.off()
  message(glue("Saved: {fname}"))

  tidy(lm_fit) %>% mutate(Catchment = catchment_label, Response = response_label)
}

# Pull selected variable names per group
sel_kn_micro <- selected_all %>%
  filter(dataset == "Kinzig k_micro") %>% pull(variable)
sel_kn_total <- selected_all %>%
  filter(dataset == "Kinzig k_total") %>% pull(variable)
sel_by_micro <- selected_all %>%
  filter(dataset == "Boye k_micro")   %>% pull(variable)
sel_by_total <- selected_all %>%
  filter(dataset == "Boye k_total")   %>% pull(variable)

# Fit and plot – Kinzig
res_kn_micro <- fit_and_plot_marginal(
  site_kinzig, "k_micro", sel_kn_micro,
  c("stream", "k_total"), "Kinzig", "k_micro",
  "Microbial k (per degree day)")

res_kn_total <- fit_and_plot_marginal(
  site_kinzig, "k_total", sel_kn_total,
  c("stream", "k_micro"), "Kinzig", "k_total",
  "Total k (per degree day)")

# Fit and plot – Boye
res_by_micro <- fit_and_plot_marginal(
  site_boye, "k_micro", sel_by_micro,
  c("stream", "k_total"), "Boye", "k_micro",
  "Microbial k (per degree day)")

res_by_total <- fit_and_plot_marginal(
  site_boye, "k_total", sel_by_total,
  c("stream", "k_micro"), "Boye", "k_total",
  "Total k (per degree day)")

# Save OLS coefficient tables
marginal_models <- bind_rows(
  res_kn_micro, res_kn_total,
  res_by_micro, res_by_total
)

if (nrow(marginal_models) > 0) {
  write.csv(marginal_models, "Marginal_effects_OLS_coefficients.csv",
            row.names = FALSE)
  message("Saved: Marginal_effects_OLS_coefficients.csv")
}


# =============================================================================
# 18  Publication-quality figures: Pattern and Mechanism outputs
# =============================================================================
#
# Two figure/table pairs designed for a master's thesis and journal submission,
# following the analytical style of Schreiner et al. (2023):
#
# Output 1 – PATTERN  (Section 18.1-18.2)
#   Question: Do decomposition rates differ between environments?
#   Best visualisation: violin + embedded boxplot + jittered raw points +
#   significance annotation. Violins reveal distribution shape that plain
#   boxplots hide. Significance brackets (2-group) and compact letter display
#   (3-group, Tukey-corrected) communicate statistical outcomes without
#   cluttering axis space.
#
# Output 2 – MECHANISM  (Section 18.3-18.4)
#   Question: What environmental factors explain these differences?
#   Best visualisation: ggpredict partial regression panels, one per predictor,
#   with stream-level observed means overlaid. This is superior to raw scatter
#   plots because it shows the marginal effect of each driver while accounting
#   for all other predictors and the random stream structure.
# =============================================================================


# ----- Shared publication theme (used by both outputs) -----------------------
#
# theme_classic removes the grey background and grid lines. We increase base
# font size, bold axis titles, and use black axis text – all standard
# requirements for journal submission (e.g., Freshwater Biology, Oecologia).

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.text        = element_text(colour = "black",   size = 11),
    axis.title       = element_text(colour = "black",   size = 12, face = "bold"),
    strip.background = element_rect(fill   = "grey95",  colour = "grey60"),
    strip.text       = element_text(face   = "bold",    size = 11),
    legend.position  = "none",
    plot.tag         = element_text(face   = "bold",    size = 14)
  )

# Colour palettes chosen for colour-blind safety (ColorBrewer diverging)
land_pal <- c("rural"  = "#2C7BB6", "urban"   = "#D7191C")
rest_pal <- c("reference"         = "#636363",
              "recently_restored" = "#3182BD",
              "long_restored"     = "#31A354")

# Helper: convert p-value to asterisk string
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}


# =============================================================================
# 18.1  Output 1 – PATTERN figure
#
# Figure title (thesis / paper):
# "Figure X. Leaf litter decomposition rates across land-use gradients and
# restoration stages. (A-B) Kinzig catchment: microbial (fine-mesh, k_micro)
# and total (coarse-mesh, k_total) decomposition rates in rural vs urban
# streams. (C-D) Boye catchment: k_micro and k_total across reference,
# recently restored (~1 yr) and long-restored (~20 yr) sites. Violins show
# kernel-density distributions; boxes show median and interquartile range;
# filled diamonds show group means; dots show individual leaf bags. Brackets
# (A-B) indicate pairwise significance from linear mixed-effects models
# (Satterthwaite p; *** < 0.001, ** < 0.01, * < 0.05, ns >= 0.05). Letters
# (C-D) indicate Tukey-corrected post-hoc groups (emmeans). k units:
# degree-day^-1."
# =============================================================================

## 18.1a Violin panel helper ---------------------------------------------------
#
# One reusable function builds a single panel: violin (distribution shape) +
# narrow boxplot (median/IQR) + jitter (raw observations) + filled diamond
# (group mean). The legend is suppressed because fill colour mirrors the x axis.

violin_panel <- function(df, x_var, y_var, palette, x_label, y_label, tag) {
  ggplot(df,
         aes(x     = .data[[x_var]],
             y     = .data[[y_var]],
             fill  = .data[[x_var]],
             colour = .data[[x_var]])) +
    # --- Distribution shape
    geom_violin(trim    = TRUE,
                alpha   = 0.32,
                colour  = NA,
                width   = 0.75) +
    # --- IQR box (white fill keeps it readable inside the violin)
    geom_boxplot(width         = 0.10,
                 outlier.shape = NA,
                 fill          = "white",
                 colour        = "grey25",
                 linewidth     = 0.55) +
    # --- Raw data points per bag
    geom_jitter(width       = 0.07,
                size        = 1.3,
                alpha       = 0.55,
                shape       = 16,
                show.legend = FALSE) +
    # --- Group mean as filled diamond
    stat_summary(fun         = mean,
                 geom        = "point",
                 shape       = 18,
                 size        = 4.5,
                 colour      = "grey10",
                 show.legend = FALSE) +
    scale_fill_manual(values  = palette) +
    scale_colour_manual(values = palette) +
    labs(x = x_label, y = y_label, tag = tag) +
    theme_pub
}


## 18.1b Significance annotations – Kinzig (2-group: significance bracket) ----
#
# For a binary predictor (rural vs urban) the geom_signif bracket with the
# mixed-model p-value is cleaner than CLD letters.

p_kn_micro_val <- tidy_kinzig_micro %>%
  filter(term == "LandUseurban") %>% pull(p.value)
p_kn_total_val <- tidy_kinzig_total %>%
  filter(term == "LandUseurban") %>% pull(p.value)


## 18.1c Compact letter display – Boye (3-group: Tukey post-hoc letters) ------
#
# emmeans computes estimated marginal means at each restoration level;
# multcomp::cld() converts Tukey contrasts into letter groups. Streams sharing
# a letter are not significantly different (p > 0.05).

get_cld_letters <- function(model, term_var) {
  # emmeans 2.x removed the cld.emmGrid method.  We compute Tukey-adjusted
  # pairwise p-values directly and feed them to multcompView::multcompLetters()
  # to obtain the compact letter display.
  em   <- emmeans::emmeans(model, as.formula(paste("~", term_var)))
  pw   <- as.data.frame(pairs(em, adjust = "tukey"))
  # multcompLetters expects a named numeric vector; names must be "A-B" format
  p_vec <- setNames(pw$p.value,
                    gsub(" ", "", gsub(" - ", "-", as.character(pw$contrast))))
  let   <- multcompView::multcompLetters(p_vec, threshold = 0.05)$Letters
  let   # named character vector: names = group levels, values = letters
}

cld_boye_micro <- get_cld_letters(model_boye_micro, "RestStatus")
cld_boye_total <- get_cld_letters(model_boye_total, "RestStatus")

# Helper: overlay CLD letters above each violin
add_cld_text <- function(p, df, y_var, letters_vec, x_var) {
  y_max  <- max(df[[y_var]], na.rm = TRUE)
  ldf    <- tibble(
    grp   = names(letters_vec),
    label = unname(letters_vec),
    y     = y_max * 1.13
  ) %>%
    mutate(grp = factor(grp, levels = levels(df[[x_var]])))
  p + geom_text(data        = ldf,
                aes(x = grp, y = y, label = label),
                inherit.aes = FALSE,
                size        = 4.5,
                fontface    = "bold",
                colour      = "grey20")
}


## 18.1d Build the four panels -------------------------------------------------

rest_labels <- c("reference"         = "Reference",
                 "recently_restored" = "Recent\n(~1 yr)",
                 "long_restored"     = "Long\n(~20 yr)")

# Panel A – Kinzig microbial
p18_A <- violin_panel(
  kinzig_micro, "LandUse", "k", land_pal,
  "Land use",
  expression(italic(k)[micro]~(degree*C^{-1}~day^{-1})),
  "A"
) +
  geom_signif(
    comparisons  = list(c("rural", "urban")),
    annotations  = sig_stars(p_kn_micro_val),
    y_position   = max(kinzig_micro$k, na.rm = TRUE) * 1.17,
    tip_length   = 0.015,
    textsize     = 5,
    vjust        = 0.3
  )

# Panel B – Kinzig total
p18_B <- violin_panel(
  kinzig_total, "LandUse", "k", land_pal,
  "Land use",
  expression(italic(k)[total]~(degree*C^{-1}~day^{-1})),
  "B"
) +
  geom_signif(
    comparisons  = list(c("rural", "urban")),
    annotations  = sig_stars(p_kn_total_val),
    y_position   = max(kinzig_total$k, na.rm = TRUE) * 1.17,
    tip_length   = 0.015,
    textsize     = 5,
    vjust        = 0.3
  )

# Panel C – Boye microbial (3-group + CLD letters)
p18_C <- violin_panel(
  boye_micro, "RestStatus", "k", rest_pal,
  "Restoration status",
  expression(italic(k)[micro]~(degree*C^{-1}~day^{-1})),
  "C"
) +
  scale_x_discrete(labels = rest_labels)
p18_C <- add_cld_text(p18_C, boye_micro, "k", cld_boye_micro, "RestStatus")

# Panel D – Boye total (3-group + CLD letters)
p18_D <- violin_panel(
  boye_total, "RestStatus", "k", rest_pal,
  "Restoration status",
  expression(italic(k)[total]~(degree*C^{-1}~day^{-1})),
  "D"
) +
  scale_x_discrete(labels = rest_labels)
p18_D <- add_cld_text(p18_D, boye_total, "k", cld_boye_total, "RestStatus")


## 18.1e Assemble and export Figure 1 (Pattern) --------------------------------
#
# Row labels are drawn as thin annotation rows using ggdraw()+draw_label().
# rel_heights gives the label rows a small fraction of the total height.

row_kn <- ggdraw() +
  draw_label("Kinzig catchment – land use effect",
             fontface = "bold", size = 10.5, colour = "grey30",
             x = 0.5, hjust = 0.5)
row_by <- ggdraw() +
  draw_label("Boye catchment – restoration status effect",
             fontface = "bold", size = 10.5, colour = "grey30",
             x = 0.5, hjust = 0.5)

fig_pattern <- plot_grid(
  row_kn,
  plot_grid(p18_A, p18_B, nrow = 1, align = "h"),
  row_by,
  plot_grid(p18_C, p18_D, nrow = 1, align = "h"),
  ncol        = 1,
  rel_heights = c(0.06, 1, 0.06, 1)
)

png("Figure1_Pattern_decomposition_rates.png",
    width = 8.5, height = 9.5, units = "in", res = 600)
print(fig_pattern)
dev.off()
message("Saved: Figure1_Pattern_decomposition_rates.png")


# =============================================================================
# 18.2  Output 1 – PATTERN model summary table
#
# Table title (thesis / paper):
# "Table X. Linear mixed-effects model results for the effect of land use
# (Kinzig catchment) and restoration status (Boye catchment) on leaf litter
# decomposition rates (k, degree-day^-1). Stream was included as a random
# intercept. Estimates represent the contrast with the reference level (rural
# or reference stream). SE = standard error; t = Satterthwaite t-statistic;
# Sig = *** p<0.001, ** p<0.01, * p<0.05, ns p>=0.05."
# =============================================================================

make_pattern_table <- function(model, catchment, endpoint, ref_label) {
  tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Catchment = catchment,
      Endpoint  = endpoint,
      Contrast  = dplyr::recode(term,
        "LandUseurban"                = glue("Urban vs {ref_label}"),
        "RestStatusrecently_restored" = glue("Recently restored vs {ref_label}"),
        "RestStatuslong_restored"     = glue("Long restored vs {ref_label}")
      ),
      Estimate = round(estimate,  6),
      SE       = round(std.error, 6),
      t        = round(statistic, 3),
      Sig      = sig_stars(p.value)
    ) %>%
    select(Catchment, Endpoint, Contrast, Estimate, SE, t, Sig)
}

table_pattern <- bind_rows(
  make_pattern_table(model_kinzig_micro, "Kinzig", "k_micro", "rural"),
  make_pattern_table(model_kinzig_total, "Kinzig", "k_total", "rural"),
  make_pattern_table(model_boye_micro,   "Boye",   "k_micro", "reference"),
  make_pattern_table(model_boye_total,   "Boye",   "k_total", "reference")
)

print(table_pattern)
write.csv(table_pattern, "Table1_Pattern_model_results.csv", row.names = FALSE)
message("Saved: Table1_Pattern_model_results.csv")


# =============================================================================
# 18.3  Output 2 – MECHANISM figure
#
# Figure title (thesis / paper):
# "Figure X. Marginal effects of environmental drivers on leaf litter
# decomposition rate from linear mixed-effects models (stream as random
# intercept). Each panel shows the predicted relationship between one
# predictor and k while holding all other predictors at their mean
# (ggpredict, ggeffects package). Shaded bands = 95% confidence intervals.
# Filled circles = stream-level mean observed k. Dashed line = overall mean k.
# Continuous predictors were standardised (z-score; x-axis in SD units).
# (A) Agricultural land use %, both catchments (H1).
# (B) Conductivity, both catchments (H1).
# (C) Sum of metal toxic units, Boye (H2).
# (D) Nitrate concentration, Kinzig (H2).
# (E) CPOM streambed cover, Kinzig (H4).
# (F) Fine sediment class (< vs >= 25% cover), Boye (H3)."
# =============================================================================

## 18.3a Partial effect panel helper -------------------------------------------
#
# ggpredict() returns a predictions dataframe with columns x, predicted,
# conf.low, conf.high. We overlay stream-level observed means (one dot per
# stream) to show how well the model tracks actual observations.
# A dashed horizontal line at the grand mean provides a visual zero-effect
# reference.

partial_panel <- function(model, pred_sc, df_obs,
                           x_label, y_label, tag, dot_colour) {
  eff    <- ggpredict(model, terms = glue("{pred_sc} [all]"))
  k_mean <- mean(df_obs$k, na.rm = TRUE)

  # Stream-level means of both response and predictor
  obs <- df_obs %>%
    group_by(stream) %>%
    summarise(
      k_site = mean(k,                      na.rm = TRUE),
      x_site = mean(.data[[pred_sc]],       na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(k_site), is.finite(x_site))

  ggplot() +
    # Confidence ribbon
    geom_ribbon(data = as.data.frame(eff),
                aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "grey70", alpha = 0.50,
                inherit.aes = FALSE) +
    # Model regression line
    geom_line(data = as.data.frame(eff),
              aes(x = x, y = predicted),
              colour = "grey20", linewidth = 1.0,
              inherit.aes = FALSE) +
    # Grand mean reference line
    geom_hline(yintercept = k_mean,
               linetype   = "dashed",
               colour     = "grey55",
               linewidth  = 0.45) +
    # Observed stream means
    geom_point(data = obs,
               aes(x = x_site, y = k_site),
               colour = dot_colour, size = 3.0,
               alpha  = 0.85, shape = 16,
               inherit.aes = FALSE) +
    labs(x = x_label, y = y_label, tag = tag) +
    theme_pub
}


## 18.3b Build five continuous panels (A-E) and one categorical panel (F) ------

# Panel A – Agricultural land use (H1, both catchments)
p18_mA <- partial_panel(
  model_H1, "pct_agri_sc", all_micro_H1,
  "Agricultural land use (%) [SD]",
  expression(italic(k)[micro]~(degree*C^{-1}~day^{-1})),
  "A", "#9E0142"
)

# Panel B – Conductivity (H1, both catchments)
p18_mB <- partial_panel(
  model_H1, "conductivity_sc", all_micro_H1,
  expression("Conductivity ("*mu*S~cm^{-1}*") [SD]"),
  expression(italic(k)[micro]~(degree*C^{-1}~day^{-1})),
  "B", "#5E4FA2"
)

# Panel C – Metal toxic units (H2, Boye)
p18_mC <- partial_panel(
  model_H2_boye, "sumTU_sc", boye_micro_H2,
  "Sum metal toxic units [SD]",
  expression(italic(k)[micro]~(degree*C^{-1}~day^{-1})),
  "C", "#D53E4F"
)

# Panel D – Nitrate (H2, Kinzig)
p18_mD <- partial_panel(
  model_H2_kinzig, "no3_sc", kinzig_micro_H2,
  expression(NO[3]^{"-"}~"concentration [SD]"),
  expression(italic(k)[micro]~(degree*C^{-1}~day^{-1})),
  "D", "#3288BD"
)

# Panel E – CPOM (H4, Kinzig)
p18_mE <- partial_panel(
  model_H4, "cpom_sc", kinzig_total_H4,
  "CPOM streambed cover (%) [SD]",
  expression(italic(k)[total]~(degree*C^{-1}~day^{-1})),
  "E", "#66C2A5"
)

# Panel F – Fine sediment class (H3, Boye)
# This is a categorical predictor so we use a violin + boxplot instead of a
# regression line. The significance bracket comes from the H3 model.

p_H3_pval <- tidy(model_H3, effects = "fixed") %>%
  filter(term == "fine_classhigh_fine") %>% pull(p.value)

p18_mF <- ggplot(boye_total_H3,
                 aes(x      = fine_class,
                     y      = k,
                     fill   = fine_class,
                     colour = fine_class)) +
  geom_violin(trim    = TRUE,  alpha   = 0.32,
              colour  = NA,    width   = 0.65) +
  geom_boxplot(width         = 0.10,
               outlier.shape = NA,
               fill          = "white",
               colour        = "grey25",
               linewidth     = 0.55) +
  geom_jitter(width = 0.07, size = 1.3,
              alpha = 0.55, shape = 16,
              show.legend = FALSE) +
  stat_summary(fun    = mean, geom = "point",
               shape  = 18,  size = 4.5,
               colour = "grey10", show.legend = FALSE) +
  geom_signif(
    comparisons  = list(c("low_fine", "high_fine")),
    annotations  = sig_stars(p_H3_pval),
    y_position   = max(boye_total_H3$k, na.rm = TRUE) * 1.17,
    tip_length   = 0.015,
    textsize     = 5,
    vjust        = 0.3
  ) +
  scale_fill_manual(values   = c("low_fine"  = "#74C476",
                                  "high_fine" = "#A1D99B")) +
  scale_colour_manual(values = c("low_fine"  = "#238B45",
                                  "high_fine" = "#006D2C")) +
  scale_x_discrete(labels = c("low_fine"  = "Low fine\n(< 25%)",
                               "high_fine" = "High fine\n(>= 25%)")) +
  labs(x   = "Fine sediment cover (Boye)",
       y   = expression(italic(k)[total]~(degree*C^{-1}~day^{-1})),
       tag = "F") +
  theme_pub


## 18.3c Assemble and export Figure 2 (Mechanism) ------------------------------

fig_mech <- plot_grid(
  p18_mA, p18_mB,
  p18_mC, p18_mD,
  p18_mE, p18_mF,
  ncol  = 2,
  nrow  = 3,
  align = "hv"
)

png("Figure2_Mechanism_environmental_drivers.png",
    width = 9, height = 12, units = "in", res = 600)
print(fig_mech)
dev.off()
message("Saved: Figure2_Mechanism_environmental_drivers.png")


# =============================================================================
# 18.4  Output 2 – MECHANISM model summary table
#
# Table title (thesis / paper):
# "Table X. Fixed-effect estimates from linear mixed-effects models testing
# the effect of environmental drivers on leaf litter decomposition rate
# (k, degree-day^-1). All continuous predictors were standardised (z-score)
# prior to modelling; estimates represent the change in k per one SD increase
# in the predictor. Positive estimates indicate faster decomposition with
# increasing predictor values. SE = standard error; t = Satterthwaite
# t-statistic; Sig = *** p<0.001, ** p<0.01, * p<0.05, ns p>=0.05;
# n_bags = leaf bags; n_sites = streams (random intercepts)."
# =============================================================================

make_mech_row <- function(model, df, predictor_map, hypothesis,
                           catchment, endpoint) {
  n_bags  <- nrow(df)
  n_sites <- n_distinct(df$stream)

  tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Hypothesis = hypothesis,
      Catchment  = catchment,
      Endpoint   = endpoint,
      Predictor  = dplyr::recode(term, !!!predictor_map),
      Estimate   = round(estimate,  6),
      SE         = round(std.error, 6),
      t          = round(statistic, 3),
      p          = ifelse(p.value < 0.001,
                          "< 0.001",
                          as.character(round(p.value, 3))),
      Sig        = sig_stars(p.value),
      n_bags     = n_bags,
      n_sites    = n_sites
    ) %>%
    select(Hypothesis, Catchment, Endpoint, Predictor,
           Estimate, SE, t, p, Sig, n_bags, n_sites)
}

table_mech <- bind_rows(

  # H1 – agriculture + conductivity (both catchments)
  make_mech_row(
    model_H1, all_micro_H1,
    c(pct_agri_sc      = "Agricultural land use (%)",
      conductivity_sc  = "Conductivity (uS/cm)"),
    "H1", "Both", "k_micro"
  ),

  # H2 – metal toxic units (Boye)
  make_mech_row(
    model_H2_boye, boye_micro_H2,
    c(sumTU_sc = "Sum metal toxic units (Cd+Cu+Fe+Ni+Pb+Zn)"),
    "H2", "Boye", "k_micro"
  ),

  # H2 – nutrients (Kinzig)
  make_mech_row(
    model_H2_kinzig, kinzig_micro_H2,
    c(no3_sc = "Nitrate (NO3, mg/L)",
      po4_sc = "Ortho-phosphate (mg/L)"),
    "H2", "Kinzig", "k_micro"
  ),

  # H3 – fine sediment class (Boye)
  make_mech_row(
    model_H3, boye_total_H3,
    c(fine_classhigh_fine = "High fine sediment (>= 25% vs < 25%)"),
    "H3", "Boye", "k_total"
  ),

  # H4 – CPOM (Kinzig)
  make_mech_row(
    model_H4, kinzig_total_H4,
    c(cpom_sc = "CPOM streambed cover (%)"),
    "H4", "Kinzig", "k_total"
  )
)

print(table_mech)
write.csv(table_mech, "Table2_Mechanism_environmental_drivers.csv",
          row.names = FALSE)
message("Saved: Table2_Mechanism_environmental_drivers.csv")


# =============================================================================
# 19  Appendix: Conceptual ecological diagrams
# =============================================================================
#
# Two publication-ready schematic figures for the thesis Appendix / Discussion.
# These are NOT raw-data plots; they represent theoretical expectations
# supported by the study's hypotheses (H1-H4) and ecological literature, and
# serve to help the reader interpret the empirical results.
#
# Figure AppC1 Panel A – Urbanisation stress gradient (Kinzig concept)
#   Shows how the three decomposition endpoints (k_micro, k_macro, k_total)
#   respond along a gradient of increasing urban / impervious cover. The
#   k_micro hump at low-moderate urbanisation reflects nutrient enrichment
#   stimulating aquatic hyphomycetes (H2); the subsequent decline reflects
#   metal and organic toxicant suppression of microbial activity (H1, H2).
#   The monotonic decline in k_macro reflects progressive shredder exclusion
#   by fine sediment accumulation (H3) and habitat simplification.
#
# Figure AppC1 Panel B – Restoration recovery trajectory (Boye concept)
#   Shows how k_micro, k_macro, and k_total recover after stream restoration.
#   k_micro recovers rapidly (microbes are r-strategists, recolonise from
#   the water column within months). k_macro follows a slower sigmoidal
#   trajectory reflecting the time required for instream habitat complexity
#   to develop and support shredding macroinvertebrate communities (H3;
#   Lepori et al. 2005). Dashed horizontal lines mark the reference state
#   (BOYkl, undisturbed) as the long-term asymptotic recovery target.
#
# Both figures use the same colour palette and line types as the data figures
# (Sections 18.1-18.3) for consistency across the thesis.
# =============================================================================


## 19.1 Shared endpoint colour / linetype key ---------------------------------
#
# These match the colours used in Section 18 so the reader can cross-reference
# between the conceptual and empirical figures without re-learning a legend.

ep_col <- c("k_total" = "#252525",
            "k_macro" = "#2166AC",
            "k_micro" = "#D6604D")

ep_lty <- c("k_total" = "solid",
            "k_macro" = "dashed",
            "k_micro" = "dotdash")

ep_lbl <- c(
  "k_total" = expression(italic(k)[total]~"(coarse mesh – total)"),
  "k_macro" = expression(italic(k)[macro]~"(macroinvert. contribution)"),
  "k_micro" = expression(italic(k)[micro]~"(fine mesh – microbial)")
)

y_axis_label <- expression(italic(k)~"("*degree*C^{-1}~day^{-1}*")")

theme_concept <- theme_classic(base_size = 11) +
  theme(
    axis.text        = element_text(colour = "black", size = 10),
    axis.title       = element_text(colour = "black", size = 11,
                                    face = "bold"),
    legend.background = element_rect(fill    = "white",
                                     colour  = "grey75",
                                     linewidth = 0.35),
    legend.key.width  = unit(1.3, "cm"),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.5),
    plot.tag          = element_text(face = "bold", size = 14),
    plot.background   = element_rect(fill = "white", colour = NA)
  )


## 19.2 Panel A – Urbanisation stress gradient (Kinzig) -----------------------
#
# We generate 400 evenly spaced x-values (0-100% urban cover) and evaluate
# three parametric curves representing the expected response of each endpoint.
#
# k_micro: Gaussian hump centred at ~18% urban (nutrient stimulation), with
#   a small linear decline at high urban cover (metal and toxicant inhibition).
# k_macro: Negative exponential decay (shredder loss driven by fine sediment
#   and habitat simplification); residual asymptote of 0.00015 represents
#   basal fungal contribution via fine-mesh process.
# k_total: Sum of the above two trajectories.
#
# Parameter choices are set to produce biologically plausible magnitudes
# (k_micro ~ 0.0018-0.0026, k_macro ~ 0.0003-0.0032) based on published
# decomposition studies in temperate streams (Graça et al. 2005;
# Woodward et al. 2012).

x_u <- seq(0, 100, length.out = 400)

k_mi_u <- 0.00200 +
  0.00045 * exp(-((x_u - 18)^2) / (2 * 14^2)) -
  0.0000022 * x_u

k_ma_u <- 0.00300 * exp(-0.033 * x_u) + 0.00015

k_to_u <- k_mi_u + k_ma_u

df_u <- bind_rows(
  tibble(x = x_u, k = k_to_u, Endpoint = "k_total"),
  tibble(x = x_u, k = k_ma_u, Endpoint = "k_macro"),
  tibble(x = x_u, k = k_mi_u, Endpoint = "k_micro")
) %>%
  mutate(Endpoint = factor(Endpoint,
                            levels = c("k_total", "k_macro", "k_micro")))

y_max_u <- max(k_to_u) * 1.30   # headroom for zone labels

# Helper to find k at a specific x value
kval <- function(kvec, x_target) kvec[which.min(abs(x_u - x_target))]

pA <- ggplot(df_u, aes(x = x, y = k,
                        colour   = Endpoint,
                        linetype = Endpoint)) +
  # ---- Zone background rectangles (plotted first = behind everything)
  annotate("rect", xmin =  0, xmax =  33,
           ymin = 0, ymax = y_max_u, fill = "#D5E8D4", alpha = 0.55) +
  annotate("rect", xmin = 33, xmax =  67,
           ymin = 0, ymax = y_max_u, fill = "#FFF2CC", alpha = 0.55) +
  annotate("rect", xmin = 67, xmax = 100,
           ymin = 0, ymax = y_max_u, fill = "#F8CECC", alpha = 0.55) +
  # ---- Zone labels (top of panel)
  annotate("text", x =  16.5, y = y_max_u * 0.955,
           label = "Rural\nreference",  size = 3.2,
           colour = "grey38", fontface = "plain") +
  annotate("text", x =  50,   y = y_max_u * 0.955,
           label = "Transitional",      size = 3.2,
           colour = "grey38", fontface = "plain") +
  annotate("text", x =  83.5, y = y_max_u * 0.955,
           label = "Urban\nimpacted",   size = 3.2,
           colour = "grey38", fontface = "plain") +
  # ---- Main trajectory curves
  geom_line(linewidth = 1.15) +
  # ---- Mechanism annotation 1: Nutrient enrichment (hump on k_micro)
  annotate("text",    x = 18, y = kval(k_mi_u, 18) + 0.00068,
           label = "Nutrient\nenrichment",
           size = 2.95, colour = "#31A354", fontface = "italic") +
  annotate("segment", x = 18, xend = 18,
           y = kval(k_mi_u, 18) + 0.00055,
           yend = kval(k_mi_u, 18) + 0.00007,
           colour = "#31A354", linewidth = 0.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  # ---- Mechanism annotation 2: Metal stress (declining k_micro)
  annotate("text",    x = 73, y = kval(k_mi_u, 73) + 0.00060,
           label = "Metal stress\n& toxicants",
           size = 2.95, colour = "#CB181D", fontface = "italic") +
  annotate("segment", x = 73, xend = 73,
           y = kval(k_mi_u, 73) + 0.00045,
           yend = kval(k_mi_u, 73) + 0.00007,
           colour = "#CB181D", linewidth = 0.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  # ---- Mechanism annotation 3: Shredder limitation (declining k_macro)
  annotate("text",    x = 83, y = kval(k_ma_u, 83) + 0.00110,
           label = "Shredder\nlimitation",
           size = 2.95, colour = "#2166AC", fontface = "italic") +
  annotate("segment", x = 83, xend = 83,
           y = kval(k_ma_u, 83) + 0.00095,
           yend = kval(k_ma_u, 83) + 0.00007,
           colour = "#2166AC", linewidth = 0.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  # ---- Scales and labels
  scale_colour_manual(values = ep_col, labels = ep_lbl) +
  scale_linetype_manual(values = ep_lty, labels = ep_lbl) +
  scale_x_continuous(breaks = seq(0, 100, 20),
                     labels = paste0(seq(0, 100, 20), "%")) +
  scale_y_continuous(
    limits = c(0, y_max_u),
    labels = function(v) formatC(v, format = "f", digits = 4)
  ) +
  labs(
    x   = "Urban / impervious land cover (%)",
    y   = y_axis_label,
    tag = "A"
  ) +
  theme_concept +
  theme(legend.position = c(0.74, 0.87))


## 19.3 Panel B – Restoration recovery trajectory (Boye) ----------------------
#
# X-axis: continuous time since restoration (0-25 yr). Vertical dotted lines
# mark the two study time points: ~1 yr (recently restored Boye sites) and
# ~20 yr (long-restored sites). Dashed horizontal lines mark the reference
# state values (BOYkl, undisturbed) as the long-term asymptotic target.
#
# k_micro: exponential approach to reference (fast recovery; fungi and bacteria
#   recolonise from the water column within months, so the trajectory is steep
#   early and flattens near the reference asymptote).
# k_macro: sigmoidal recovery (S-shaped; slow initial establishment of
#   shredders limited by instream wood and hydraulic complexity, then faster
#   growth as habitat matures, then deceleration near the reference asymptote).

t_r <- seq(0, 25, length.out = 400)

# Reference (target) asymptotic values for each endpoint
ref_mi <- 0.00220
ref_ma <- 0.00300
ref_to <- ref_mi + ref_ma

k_mi_r <- ref_mi * (1 - 0.28 * exp(-0.50 * t_r))
k_ma_r <- ref_ma / (1 + exp(-0.28 * (t_r - 9)))
k_to_r <- k_mi_r + k_ma_r

df_r <- bind_rows(
  tibble(t = t_r, k = k_to_r, Endpoint = "k_total"),
  tibble(t = t_r, k = k_ma_r, Endpoint = "k_macro"),
  tibble(t = t_r, k = k_mi_r, Endpoint = "k_micro")
) %>%
  mutate(Endpoint = factor(Endpoint,
                            levels = c("k_total", "k_macro", "k_micro")))

y_max_r <- ref_to * 1.30
kvalr   <- function(kvec, t_target) kvec[which.min(abs(t_r - t_target))]

pB <- ggplot(df_r, aes(x = t, y = k,
                        colour   = Endpoint,
                        linetype = Endpoint)) +
  # ---- Zone backgrounds
  annotate("rect", xmin =  0, xmax =  4,
           ymin = 0, ymax = y_max_r, fill = "#F8CECC", alpha = 0.55) +
  annotate("rect", xmin =  4, xmax = 14,
           ymin = 0, ymax = y_max_r, fill = "#FFF2CC", alpha = 0.55) +
  annotate("rect", xmin = 14, xmax = 25,
           ymin = 0, ymax = y_max_r, fill = "#D5E8D4", alpha = 0.55) +
  # ---- Zone labels
  annotate("text", x =  2,    y = y_max_r * 0.955,
           label = "Early\nrecovery",         size = 3.2, colour = "grey38") +
  annotate("text", x =  9,    y = y_max_r * 0.955,
           label = "Progressive\nrecovery",   size = 3.2, colour = "grey38") +
  annotate("text", x = 19.5,  y = y_max_r * 0.955,
           label = "Advanced\nrecovery",      size = 3.2, colour = "grey38") +
  # ---- Study time-point markers
  annotate("segment", x = 1, xend = 1,
           y = 0, yend = y_max_r * 0.80,
           linetype = "dotted", colour = "grey50", linewidth = 0.70) +
  annotate("text", x = 1, y = y_max_r * 0.84,
           label = "~1 yr\n(BOY recent)", size = 2.65,
           colour = "grey45", fontface = "plain") +
  annotate("segment", x = 20, xend = 20,
           y = 0, yend = y_max_r * 0.80,
           linetype = "dotted", colour = "grey50", linewidth = 0.70) +
  annotate("text", x = 20, y = y_max_r * 0.84,
           label = "~20 yr\n(BOY long)", size = 2.65,
           colour = "grey45", fontface = "plain") +
  # ---- Reference asymptote dashed lines (drawn BEFORE curves so curves sit on top)
  geom_hline(yintercept = ref_to, linetype = "longdash",
             colour = ep_col["k_total"], linewidth = 0.50, alpha = 0.50) +
  geom_hline(yintercept = ref_ma, linetype = "longdash",
             colour = ep_col["k_macro"], linewidth = 0.50, alpha = 0.50) +
  geom_hline(yintercept = ref_mi, linetype = "longdash",
             colour = ep_col["k_micro"], linewidth = 0.50, alpha = 0.50) +
  annotate("text", x = 24.2, y = ref_to + 0.000060,
           label = "Reference\nstate (BOYkl)",
           size = 2.65, colour = "grey42", fontface = "italic", hjust = 1) +
  # ---- Main trajectory curves
  geom_line(linewidth = 1.15) +
  # ---- Mechanism annotation 1: rapid microbial recovery
  annotate("text",    x = 3.0, y = kvalr(k_mi_r, 1) + 0.00075,
           label = "Rapid\nmicrobial recovery",
           size = 2.95, colour = "#D6604D", fontface = "italic") +
  annotate("segment", x = 1.8, xend = 1.1,
           y = kvalr(k_mi_r, 1.8) + 0.00048,
           yend = kvalr(k_mi_r, 1.1) + 0.00010,
           colour = "#D6604D", linewidth = 0.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  # ---- Mechanism annotation 2: shredder recolonisation (sigmoid midpoint)
  annotate("text",    x = 9.5, y = kvalr(k_ma_r, 9) + 0.00085,
           label = "Shredder\nrecolonisation",
           size = 2.95, colour = "#2166AC", fontface = "italic") +
  annotate("segment", x = 9.0, xend = 9.0,
           y = kvalr(k_ma_r, 9) + 0.00068,
           yend = kvalr(k_ma_r, 9) + 0.00010,
           colour = "#2166AC", linewidth = 0.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  # ---- Mechanism annotation 3: restoration recovery (broad label)
  annotate("text", x = 17.5, y = y_max_r * 0.67,
           label = "Restoration\nrecovery",
           size = 3.1, colour = "#525252", fontface = "bold.italic") +
  # ---- Scales and labels
  scale_colour_manual(values = ep_col, labels = ep_lbl) +
  scale_linetype_manual(values = ep_lty, labels = ep_lbl) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25),
                     labels = paste0(c(0, 5, 10, 15, 20, 25), " yr")) +
  scale_y_continuous(
    limits = c(0, y_max_r),
    labels = function(v) formatC(v, format = "f", digits = 4)
  ) +
  labs(
    x   = "Time since restoration (years)",
    y   = y_axis_label,
    tag = "B"
  ) +
  theme_concept +
  theme(legend.position = c(0.22, 0.82))


## 19.4 Assemble, export, and print caption ------------------------------------

fig_AppC1 <- plot_grid(
  pA, pB,
  nrow        = 1,
  align       = "hv",
  rel_widths  = c(1, 1)
)

png("Figure_AppC1_Conceptual_diagrams.png",
    width = 13, height = 5.5, units = "in", res = 600)
print(fig_AppC1)
dev.off()
message("Saved: Figure_AppC1_Conceptual_diagrams.png")


## 19.5 Suggested thesis caption -----------------------------------------------

cat(
  "\n",
  "SUGGESTED APPENDIX CAPTION\n",
  "============================================================\n",
  "Figure AppC1. Conceptual diagrams of leaf litter decomposition responses\n",
  "along environmental stress and restoration gradients.\n\n",
  "(A) Urbanisation gradient (Kinzig catchment). Predicted trajectories for\n",
  "microbial (k_micro, fine-mesh), total (k_total, coarse-mesh), and\n",
  "macroinvertebrate-mediated (k_macro = k_total - k_micro) leaf litter\n",
  "decomposition rates along a gradient of increasing urban/impervious\n",
  "land cover. The hump in k_micro at low-moderate urbanisation reflects\n",
  "stimulation of aquatic hyphomycetes by moderate nutrient enrichment\n",
  "(Gulis & Suberkropp 2003). Subsequent decline at high urban cover\n",
  "reflects metal contamination and micropollutant toxicity reducing\n",
  "microbial decomposer activity (Duarte et al. 2008; Ferreira et al. 2016).\n",
  "The monotonic decline in k_macro reflects progressive shredder\n",
  "exclusion by fine sediment accumulation and habitat simplification\n",
  "(Lepori et al. 2005).\n\n",
  "(B) Restoration recovery trajectory (Boye catchment). Predicted\n",
  "trajectories of the three decomposition endpoints following stream\n",
  "restoration. k_micro recovers rapidly (exponential; microbes recolonise\n",
  "from the water column within weeks to months). k_macro follows a\n",
  "sigmoidal trajectory reflecting the slower re-establishment of shredding\n",
  "macroinvertebrates as instream habitat complexity develops (Pascoal &\n",
  "Cassio 2004). Dashed horizontal lines indicate the undisturbed reference\n",
  "state (BOYkl) as the long-term asymptotic recovery target. Dotted\n",
  "vertical lines mark the two study time points (~1 yr: recently restored\n",
  "sites; ~20 yr: long-restored sites). Curves are schematic and based on\n",
  "theoretical expectations from H1-H4 and published literature; see\n",
  "Figures 1-2 for empirical results.\n\n",
  "Colour and line-type coding is consistent across all thesis figures:\n",
  "black solid = k_total; blue dashed = k_macro; red dot-dash = k_micro.\n",
  "============================================================\n"
)


# =============================================================================
# 20  Extended analysis – setup (colours, theme, output path)
# =============================================================================

## 20.1 Colour palette (fixed throughout extension) ----------------------------
COL_RURAL     <- "#2E7D32"
COL_URBAN     <- "#E64A19"
COL_LOW_FINE  <- "#4C78A8"
COL_HIGH_FINE <- "#F58518"
COL_PAL_LU    <- c(rural    = COL_RURAL,    urban     = COL_URBAN)
COL_PAL_FIN   <- c(low_fine = COL_LOW_FINE, high_fine = COL_HIGH_FINE)

## 20.2 Shared publication theme -----------------------------------------------
theme_pub <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(
      strip.background       = element_blank(),
      panel.grid.major.y     = element_line(colour = "grey92", linewidth = 0.3),
      axis.line              = element_line(colour = "grey30"),
      axis.ticks             = element_line(colour = "grey30"),
      plot.title             = element_text(size = base,     face = "bold"),
      plot.subtitle          = element_text(size = base - 2, colour = "grey40"),
      legend.position        = "bottom",
      legend.key.size        = unit(0.4, "cm")
    )
}

## 20.3 Output path (same working directory as the rest of the script) ---------
# getwd() returns the data/processed directory set by setwd() above.
# All plots and tables are written there, consistent with Figures 1–4.
ext_out <- getwd()

save_fig <- function(plot_obj, fname, w = 16, h = 10, dpi = 300) {
  ggsave(file.path(ext_out, fname), plot_obj,
         width = w, height = h, units = "cm", dpi = dpi)
  message("Saved: ", fname)
}


# =============================================================================
# 21  H1 modern plots  (agricultural land use & conductivity → k_micro)
# =============================================================================
message("=== Section 21: H1 modern plots ===")

## 21.1 Style 1 – violin + boxplot + jitter + mean by agri tertile -------------
# Bin pct_agri into three equal-frequency groups for categorical violin display
all_micro_H1 <- all_micro_H1 %>%
  mutate(agri_tertile = cut(pct_agri,
                            breaks = quantile(pct_agri, probs = c(0, 1/3, 2/3, 1),
                                              na.rm = TRUE),
                            labels = c("Low agri", "Med agri", "High agri"),
                            include.lowest = TRUE))

p_H1_s1 <- ggplot(all_micro_H1,
                   aes(agri_tertile, k, fill = agri_tertile, colour = agri_tertile)) +
  geom_violin(alpha = 0.25, trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.16, alpha = 0.7, outlier.shape = NA,
               colour = "grey30", fill = "white") +
  geom_jitter(width = 0.1, size = 1.4, alpha = 0.55) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "white", colour = "grey20") +
  scale_fill_manual(values   = c("Low agri" = COL_RURAL,
                                  "Med agri" = "#81C784",
                                  "High agri" = COL_URBAN),
                    guide = "none") +
  scale_colour_manual(values = c("Low agri" = COL_RURAL,
                                  "Med agri" = "#81C784",
                                  "High agri" = COL_URBAN),
                      guide = "none") +
  labs(title = "H1 – Microbial decomposition by agricultural land-use class",
       x = NULL, y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 21.2 Style 2 – clean minimal effect plots (Nature / GCB style) --------------
p_H1_s2_agri <- ggplot(as.data.frame(eff_H1_agri), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_RURAL, alpha = 0.15) +
  geom_line(colour = COL_RURAL, linewidth = 1.2) +
  geom_rug(data = all_micro_H1, aes(x = pct_agri_sc),
           sides = "b", alpha = 0.35, colour = COL_RURAL,
           inherit.aes = FALSE) +
  labs(title    = "H1 – Agricultural land use effect on microbial decomposition",
       subtitle = "Partial effect ± 95 % CI (lmer, random intercept per stream)",
       x = "Agricultural land use (standardised)",
       y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

p_H1_s2_cond <- ggplot(as.data.frame(eff_H1_cond), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_URBAN, alpha = 0.15) +
  geom_line(colour = COL_URBAN, linewidth = 1.2) +
  geom_rug(data = all_micro_H1, aes(x = conductivity_sc),
           sides = "b", alpha = 0.35, colour = COL_URBAN,
           inherit.aes = FALSE) +
  labs(title    = "H1 – Conductivity effect on microbial decomposition",
       subtitle = "Partial effect ± 95 % CI",
       x = "Conductivity (standardised)",
       y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 21.3 Style 3 – annotated coefficient plot -----------------------------------
tidy_H1_ext <- broom.mixed::tidy(model_H1, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term_lab = recode(term,
    pct_agri_sc     = "Agricultural\nland use",
    conductivity_sc = "Conductivity"))

p_H1_s3 <- ggplot(tidy_H1_ext, aes(estimate, term_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.12, linewidth = 0.8, colour = "grey40") +
  geom_point(aes(colour = estimate > 0), size = 4) +
  scale_colour_manual(values  = c(`TRUE` = COL_URBAN, `FALSE` = COL_RURAL),
                      labels  = c("Positive effect", "Negative effect"),
                      name    = NULL) +
  labs(title    = "H1 – Fixed-effect coefficients",
       subtitle = "Estimate ± 95 % CI; standardised predictors",
       x = "Standardised coefficient", y = NULL) +
  theme_pub()

save_fig(p_H1_s1,       "Figure_H1_style1_violin.png")
save_fig(p_H1_s2_agri,  "Figure_H1_style2_agri.png")
save_fig(p_H1_s2_cond,  "Figure_H1_style2_cond.png")
save_fig(p_H1_s3,       "Figure_H1_style3_coefplot.png", h = 8)


# =============================================================================
# 22  H2 modern plots  (metals → k_micro Boye; nutrients → k_micro Kinzig)
# =============================================================================
message("=== Section 22: H2 modern plots ===")

## 22.1 Style 1 – violin + jitter by ΣTU class (Boye) -------------------------
boye_micro_H2_plot <- boye_micro_H2 %>%
  mutate(TU_group = factor(
    if_else(sumTU >= median(sumTU, na.rm = TRUE), "High ΣTU", "Low ΣTU"),
    levels = c("Low ΣTU", "High ΣTU")))

p_H2_s1 <- ggplot(boye_micro_H2_plot,
                   aes(TU_group, k, fill = TU_group, colour = TU_group)) +
  geom_violin(alpha = 0.25, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA,
               colour = "grey30", fill = "white") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "white", colour = "grey20") +
  scale_fill_manual(values   = c("Low ΣTU" = COL_LOW_FINE,
                                  "High ΣTU" = COL_HIGH_FINE),
                    guide = "none") +
  scale_colour_manual(values = c("Low ΣTU" = COL_LOW_FINE,
                                  "High ΣTU" = COL_HIGH_FINE),
                      guide = "none") +
  labs(title = "H2 – Microbial decomposition by metal toxic-unit class (Boye)",
       x = NULL, y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 22.2 Style 2 – marginal effect plots ----------------------------------------
p_H2_s2_tu <- ggplot(as.data.frame(eff_H2_boye), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_HIGH_FINE, alpha = 0.15) +
  geom_line(colour = COL_HIGH_FINE, linewidth = 1.2) +
  geom_rug(data = boye_micro_H2, aes(x = sumTU_sc),
           sides = "b", alpha = 0.35, colour = COL_HIGH_FINE,
           inherit.aes = FALSE) +
  labs(title    = "H2 – Metal toxic units → microbial decomposition (Boye)",
       subtitle = "Partial effect ± 95 % CI",
       x = expression(Sigma*"TU (standardised)"),
       y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

p_H2_s2_no3 <- ggplot(as.data.frame(eff_H2_no3), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_RURAL, alpha = 0.15) +
  geom_line(colour = COL_RURAL, linewidth = 1.2) +
  geom_rug(data = kinzig_micro_H2, aes(x = no3_sc),
           sides = "b", alpha = 0.35, colour = COL_RURAL,
           inherit.aes = FALSE) +
  labs(title    = expression("H2 – NO"[3]^{"-"}~"→ microbial decomposition (Kinzig)"),
       subtitle = "Partial effect ± 95 % CI",
       x = expression("NO"[3]^{"-"}~"(standardised)"),
       y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

p_H2_s2_po4 <- ggplot(as.data.frame(eff_H2_po4), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_URBAN, alpha = 0.15) +
  geom_line(colour = COL_URBAN, linewidth = 1.2) +
  geom_rug(data = kinzig_micro_H2, aes(x = po4_sc),
           sides = "b", alpha = 0.35, colour = COL_URBAN,
           inherit.aes = FALSE) +
  labs(title    = expression("H2 – PO"[4]^{"3-"}~"→ microbial decomposition (Kinzig)"),
       subtitle = "Partial effect ± 95 % CI",
       x = expression("PO"[4]^{"3-"}~"(standardised)"),
       y = expression(italic(k)[micro] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 22.3 Style 3 – coefficient comparison plot ----------------------------------
tidy_H2_b <- broom.mixed::tidy(model_H2_boye,   effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(dataset  = "Boye (ΣTU)",
         term_lab = "Σ Toxic Units")

tidy_H2_k <- broom.mixed::tidy(model_H2_kinzig, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(dataset  = "Kinzig (nutrients)",
         term_lab = recode(term, no3_sc = "NO₃⁻", po4_sc = "PO₄³⁻"))

p_H2_s3 <- bind_rows(tidy_H2_b, tidy_H2_k) %>%
  ggplot(aes(estimate, term_lab, colour = dataset)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.15, linewidth = 0.7,
                 position = position_dodge(0.35)) +
  geom_point(size = 3.5, position = position_dodge(0.35)) +
  scale_colour_manual(values = c("Boye (ΣTU)"       = COL_HIGH_FINE,
                                  "Kinzig (nutrients)" = COL_RURAL),
                      name = "System") +
  labs(title    = "H2 – Fixed-effect coefficients (Boye & Kinzig)",
       subtitle = "Estimate ± 95 % CI; standardised predictors",
       x = "Standardised coefficient", y = NULL) +
  theme_pub()

save_fig(p_H2_s1,       "Figure_H2_style1_violin_TU.png")
save_fig(p_H2_s2_tu,    "Figure_H2_style2_TU.png")
save_fig(p_H2_s2_no3,   "Figure_H2_style2_NO3.png")
save_fig(p_H2_s2_po4,   "Figure_H2_style2_PO4.png")
save_fig(p_H2_s3,       "Figure_H2_style3_coefplot.png", h = 8)


# =============================================================================
# 23  H3 modern plots  (fine sediment → k_total, Boye)
# =============================================================================
message("=== Section 23: H3 modern plots ===")

## 23.1 Style 1 – violin + boxplot + jitter by fine_class ----------------------
p_H3_s1 <- ggplot(boye_total_H3,
                   aes(fine_class, k, fill = fine_class, colour = fine_class)) +
  geom_violin(alpha = 0.25, trim = FALSE) +
  geom_boxplot(width = 0.16, alpha = 0.7, outlier.shape = NA,
               colour = "grey30", fill = "white") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "white", colour = "grey20") +
  scale_fill_manual(values   = COL_PAL_FIN, guide = "none") +
  scale_colour_manual(values = COL_PAL_FIN, guide = "none") +
  scale_x_discrete(labels = c(low_fine  = "Low fine cover\n(< 25 %)",
                               high_fine = "High fine cover\n(≥ 25 %)")) +
  labs(title = "H3 – Total decomposition by fine-sediment class (Boye)",
       x = NULL,
       y = expression(italic(k)[total] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 23.2 Style 2 – continuous fine_cover marginal effect ------------------------
model_H3_cont <- lmer(k ~ fine_cover + (1 | stream), data = boye_total_H3)
eff_H3_cont   <- ggpredict(model_H3_cont, terms = "fine_cover [all]")

p_H3_s2 <- ggplot(as.data.frame(eff_H3_cont), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_HIGH_FINE, alpha = 0.15) +
  geom_line(colour = COL_HIGH_FINE, linewidth = 1.2) +
  geom_rug(data = boye_total_H3, aes(x = fine_cover),
           sides = "b", alpha = 0.35, colour = COL_HIGH_FINE,
           inherit.aes = FALSE) +
  geom_vline(xintercept = 25, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 27, y = max(eff_H3_cont$conf.high) * 0.97,
           label = "25 % threshold", size = 3, hjust = 0, colour = "grey40") +
  labs(title    = "H3 – Fine sediment cover → total decomposition (Boye)",
       subtitle = "Partial effect ± 95 % CI (continuous lmer model)",
       x        = "Fine cover (FPOM + sand, %)",
       y        = expression(italic(k)[total] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 23.3 Style 3 – annotated boxplot with CLD letters --------------------------
# derive p-value for annotation; fall back gracefully if too few observations
H3_ann <- tryCatch({
  em_H3  <- emmeans::emmeans(model_H3, ~ fine_class)
  pw_H3  <- as.data.frame(pairs(em_H3, adjust = "tukey"))
  p_val  <- pw_H3$p.value[1]
  letter_df <- boye_total_H3 %>%
    group_by(fine_class) %>%
    summarise(y_pos = max(k, na.rm = TRUE) * 1.08, .groups = "drop") %>%
    mutate(letter = if_else(p_val < 0.05,
                             c(low_fine = "a", high_fine = "b")[as.character(fine_class)],
                             "a"))
  letter_df
}, error = function(e) NULL)

p_H3_s3 <- ggplot(boye_total_H3, aes(fine_class, k)) +
  geom_boxplot(aes(fill = fine_class), width = 0.4, alpha = 0.65,
               outlier.shape = 21, outlier.alpha = 0.5) +
  geom_jitter(aes(colour = fine_class), width = 0.12, size = 1.5, alpha = 0.5) +
  { if (!is.null(H3_ann))
      geom_text(data = H3_ann, aes(y = y_pos, label = letter),
                size = 5, fontface = "bold") } +
  scale_fill_manual(values   = COL_PAL_FIN, guide = "none") +
  scale_colour_manual(values = COL_PAL_FIN, guide = "none") +
  scale_x_discrete(labels = c(low_fine  = "Low fine\n(< 25 %)",
                               high_fine = "High fine\n(≥ 25 %)")) +
  annotate("text", x = 1.5, y = max(boye_total_H3$k, na.rm = TRUE) * 1.17,
           label = "Fine sediment may suppress\nmacroinvertebrate access",
           size = 3.1, colour = "grey35", hjust = 0.5) +
  labs(title = "H3 – Fine sediment suppresses total decomposition (Boye)",
       x = NULL, y = expression(italic(k)[total] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

save_fig(p_H3_s1, "Figure_H3_style1_violin.png")
save_fig(p_H3_s2, "Figure_H3_style2_continuous.png")
save_fig(p_H3_s3, "Figure_H3_style3_annotated.png")


# =============================================================================
# 24  H4 modern plots  (CPOM → k_total, Kinzig)
# =============================================================================
message("=== Section 24: H4 modern plots ===")

## 24.1 Style 1 – violin + boxplot + jitter by LandUse -------------------------
p_H4_s1 <- ggplot(kinzig_total_H4,
                   aes(LandUse, k, fill = LandUse, colour = LandUse)) +
  geom_violin(alpha = 0.25, trim = FALSE) +
  geom_boxplot(width = 0.16, alpha = 0.7, outlier.shape = NA,
               colour = "grey30", fill = "white") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "white", colour = "grey20") +
  scale_fill_manual(values   = COL_PAL_LU, guide = "none") +
  scale_colour_manual(values = COL_PAL_LU, guide = "none") +
  labs(title = "H4 – Total decomposition by land use (Kinzig)",
       x = NULL, y = expression(italic(k)[total] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 24.2 Style 2 – CPOM marginal effect -----------------------------------------
p_H4_s2 <- ggplot(as.data.frame(eff_H4), aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = COL_RURAL, alpha = 0.15) +
  geom_line(colour = COL_RURAL, linewidth = 1.2) +
  geom_rug(data = kinzig_total_H4, aes(x = cpom_sc),
           sides = "b", alpha = 0.35, colour = COL_RURAL,
           inherit.aes = FALSE) +
  labs(title    = "H4 – CPOM coverage → total decomposition (Kinzig)",
       subtitle = "Partial effect ± 95 % CI",
       x        = "CPOM coverage (standardised)",
       y        = expression(italic(k)[total] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

## 24.3 Style 3 – scatter + per-LandUse OLS trend + ecology annotation ---------
p_H4_s3 <- ggplot(kinzig_total_H4, aes(cpom_sc, k, colour = LandUse,
                                         fill = LandUse)) +
  geom_point(size = 2.2, alpha = 0.65) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 1) +
  scale_colour_manual(values = COL_PAL_LU, name = "Land use") +
  scale_fill_manual(values   = COL_PAL_LU, name = "Land use") +
  annotate("text",
           x = max(kinzig_total_H4$cpom_sc, na.rm = TRUE) * 0.85,
           y = max(kinzig_total_H4$k,       na.rm = TRUE) * 0.92,
           label = "More CPOM → more\nshredder activity",
           size = 3.1, colour = "grey35", hjust = 1) +
  labs(title    = "H4 – CPOM as allochthonous subsidy (Kinzig)",
       subtitle = "OLS trend per land-use class ± SE",
       x        = "CPOM coverage (standardised)",
       y        = expression(italic(k)[total] ~ (degree*C^{-1}~d^{-1}))) +
  theme_pub()

save_fig(p_H4_s1, "Figure_H4_style1_violin.png")
save_fig(p_H4_s2, "Figure_H4_style2_marginal.png")
save_fig(p_H4_s3, "Figure_H4_style3_scatter.png")


# =============================================================================
# 25  Multi-panel patchwork figures (H1→A+B, H2→A+B+C, H3, H4→A+B)
# =============================================================================
message("=== Section 25: Multi-panel patchwork figures ===")

panel_H1 <- (p_H1_s2_agri | p_H1_s2_cond) +
  plot_annotation(
    tag_levels = "A",
    title      = "H1 – Urban stress: agricultural land use and conductivity",
    theme      = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(panel_H1, "Figure_H1_multipanel.png", w = 30, h = 12)

panel_H2 <- (p_H2_s2_tu | p_H2_s2_no3 | p_H2_s2_po4) +
  plot_annotation(
    tag_levels = "A",
    title      = "H2 – Chemical stress: metal toxic units (Boye) & nutrients (Kinzig)",
    theme      = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(panel_H2, "Figure_H2_multipanel.png", w = 42, h = 12)

panel_H3 <- p_H3_s2 +
  plot_annotation(
    title = "H3 – Physical stress: fine sediment clogging (Boye)",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(panel_H3, "Figure_H3_panel.png", w = 16, h = 12)

panel_H4 <- (p_H4_s2 | p_H4_s3) +
  plot_annotation(
    tag_levels = "A",
    title      = "H4 – Resource subsidy: CPOM coverage → total decomposition (Kinzig)",
    theme      = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(panel_H4, "Figure_H4_multipanel.png", w = 30, h = 12)


# =============================================================================
# 26  Similarity / opposition analysis (k endpoints)
# =============================================================================
message("=== Section 26: Similarity / opposition analysis ===")

## 26.1 Helper: fit lmer per endpoint and extract slope sign -------------------
fit_endpoint_sign <- function(df, predictor_col, endpoint_label) {
  df2 <- df %>%
    select(k, stream, val = all_of(predictor_col)) %>%
    drop_na() %>%
    mutate(val_sc = scale(val)[, 1])
  if (nrow(df2) < 6) return(tibble(endpoint = endpoint_label, estimate = NA_real_))
  m <- tryCatch(
    lmer(k ~ val_sc + (1 | stream), data = df2),
    error = function(e) NULL)
  if (is.null(m)) return(tibble(endpoint = endpoint_label, estimate = NA_real_))
  tibble(endpoint = endpoint_label, estimate = fixef(m)["val_sc"])
}

## 26.2 Build sign matrix for pct_agri and conductivity across three endpoints -
# Kinzig endpoints: k_micro, k_total; macro uses k_macro column
kinzig_macro_sign <- kinzig_macro %>%
  select(stream, LandUse, k = k_macro) %>%
  left_join(master_kinzig %>%
              select(stream, pct_agri, conductivity_u_s_cm) %>% distinct(),
            by = "stream")

build_sim_row <- function(pred_col, pred_label) {
  r_micro <- fit_endpoint_sign(kinzig_micro, pred_col, "k_micro")
  r_total <- fit_endpoint_sign(kinzig_total, pred_col, "k_total")
  r_macro <- fit_endpoint_sign(kinzig_macro_sign, pred_col, "k_macro")
  signs   <- c(sign(r_micro$estimate), sign(r_total$estimate),
               sign(r_macro$estimate))
  signs   <- signs[!is.na(signs)]
  pct     <- if (length(signs) >= 2)
               sum(combn(signs, 2, FUN = function(p) p[1] == p[2])) /
               ncol(combn(signs, 2)) * 100 else NA_real_
  tibble(
    predictor       = pred_label,
    beta_micro      = r_micro$estimate,
    beta_total      = r_total$estimate,
    beta_macro      = r_macro$estimate,
    sign_micro      = sign(r_micro$estimate),
    sign_total      = sign(r_total$estimate),
    sign_macro      = sign(r_macro$estimate),
    similarity_pct  = pct,
    direction_class = case_when(
      pct == 100 ~ "Similar",
      pct ==   0 ~ "Opposing",
      TRUE       ~ "Mixed")
  )
}

sim_table <- bind_rows(
  build_sim_row("pct_agri",          "Agricultural land use"),
  build_sim_row("conductivity_u_s_cm","Conductivity")
)

write.csv(sim_table, "Table_similarity_k_endpoints.csv", row.names = FALSE)
message("Saved: Table_similarity_k_endpoints.csv")

## 26.3 Heatmap – effect direction per endpoint --------------------------------
sim_long <- sim_table %>%
  select(predictor, sign_micro, sign_total, sign_macro) %>%
  pivot_longer(-predictor, names_to = "endpoint", values_to = "sign_val") %>%
  mutate(
    endpoint = recode(endpoint,
                      sign_micro = "k[micro]",
                      sign_total = "k[total]",
                      sign_macro = "k[macro]"),
    sign_lab = case_when(sign_val >  0 ~ "+",
                         sign_val <  0 ~ "−",
                         TRUE          ~ "0"))

p_sim_heat <- ggplot(sim_long, aes(endpoint, predictor, fill = factor(sign_val))) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sign_lab), size = 6, fontface = "bold") +
  scale_fill_manual(
    values   = c(`-1` = COL_URBAN, `0` = "#BDBDBD", `1` = COL_RURAL),
    labels   = c("Negative", "Zero", "Positive"),
    name     = "Effect direction",
    na.value = "grey88") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  labs(title = "Similarity of predictor effects across k endpoints",
       x = "Decomposition endpoint", y = "Predictor") +
  theme_pub() + theme(legend.position = "right")

save_fig(p_sim_heat, "Figure_similarity_heatmap.png", h = 8)

## 26.4 Barplot – similarity % -------------------------------------------------
p_sim_bar <- sim_table %>%
  filter(!is.na(similarity_pct)) %>%
  mutate(predictor = fct_reorder(predictor, similarity_pct)) %>%
  ggplot(aes(predictor, similarity_pct, fill = direction_class)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = c(Similar  = COL_RURAL,
                                Opposing = COL_URBAN,
                                Mixed    = COL_HIGH_FINE),
                    name = "Direction class") +
  coord_flip() +
  ylim(0, 100) +
  labs(title    = "Predictor similarity across k endpoints",
       subtitle = "100 % = same effect direction in k_micro, k_total, k_macro",
       x = NULL, y = "Similarity (%)") +
  theme_pub()

save_fig(p_sim_bar, "Figure_similarity_barplot.png", h = 8)


# =============================================================================
# 27  Ecological interpretation table (Part 7)
# =============================================================================
message("=== Section 27: Ecological interpretation table ===")

eco_interp <- tribble(
  ~predictor,              ~direction_class, ~interpretation,
  "Agricultural land use", "Similar",
    "Consistent negative effect: elevated fine sediment, pesticide inputs, and reduced riparian shade collectively suppress aquatic hyphomycetes and bacteria.",
  "Conductivity",          "Mixed",
    "Context-dependent: moderate conductivity may stimulate microbial decomposers via nutrient co-variation, while extreme values in urban streams inhibit decomposition.",
  "Sigma Toxic Units",     "Similar",
    "Negative across endpoints: metal contamination in the Boye catchment exerts direct toxicity on fungal and bacterial decomposers (Duarte et al. 2008).",
  "NO3",                   "Mixed",
    "Low NO3 stimulates microbial activity (nutrient subsidy to fungi); high urban loads may reverse this via co-occurring toxicants (Gulis & Suberkropp 2003).",
  "PO4",                   "Mixed",
    "Orthophosphate effects are context-dependent; stimulatory at low background concentrations but potentially masked by co-stressors in urban Kinzig streams.",
  "Fine sediment cover",   "Similar",
    "Consistent negative effect on k_total: fine particles clog hyporheic interstices and exclude macroinvertebrate shredders (Lepori et al. 2005).",
  "CPOM coverage",         "Similar",
    "Consistent positive effect: CPOM is a resource subsidy stimulating shredder invertebrates and indirectly enhancing microbial colonisation (Graca et al. 2001)."
)

write.csv(eco_interp, "Table_ecological_interpretation.csv", row.names = FALSE)
message("Saved: Table_ecological_interpretation.csv")
# =============================================================================
# 29  DLM response diagnostics (Step 2)
# =============================================================================
message("=== Section 29: DLM response diagnostics ===")

dlm_list <- list(
  kinzig_micro = kinzig_micro_dlm$DLM,
  kinzig_total = kinzig_total_dlm$DLM,
  boye_micro   = boye_micro_dlm$DLM,
  boye_total   = boye_total_dlm$DLM
)

## 29.1 Skewness table ---------------------------------------------------------
skew_dlm <- tibble(
  endpoint        = names(dlm_list),
  n               = sapply(dlm_list, length),
  mean_dlm        = sapply(dlm_list, mean,            na.rm = TRUE),
  sd_dlm          = sapply(dlm_list, sd,              na.rm = TRUE),
  skewness        = sapply(dlm_list, e1071::skewness, na.rm = TRUE),
  log10_transform = abs(sapply(dlm_list, e1071::skewness, na.rm = TRUE)) > 2
)
write.csv(skew_dlm, "Table_DLM_skewness.csv", row.names = FALSE)

## 29.2 Histogram panel --------------------------------------------------------
hist_plots <- lapply(names(dlm_list), function(nm) {
  sk <- round(e1071::skewness(dlm_list[[nm]], na.rm = TRUE), 2)
  tibble(DLM = dlm_list[[nm]]) %>%
    ggplot(aes(DLM)) +
    geom_histogram(bins = 18, fill = COL_RURAL, colour = "white", alpha = 0.85) +
    labs(title    = nm,
         subtitle = glue("skewness = {sk}"),
         x = "DLM (% °C⁻¹ d⁻¹)", y = "Count") +
    theme_pub(9)
})

p_dlm_hist <- wrap_plots(hist_plots, ncol = 2) +
  plot_annotation(title = "DLM response distributions by endpoint")
save_fig(p_dlm_hist, "Figure_DLM_histograms.png", w = 28, h = 18)

## 29.3 Pearson correlation micro vs total (Kinzig) ----------------------------
cor_kz_dlm <- tryCatch({
  df_cor <- kinzig_micro_dlm %>%
    select(stream, replicate, DLM_micro = DLM) %>%
    inner_join(kinzig_total_dlm %>% select(stream, replicate, DLM_total = DLM),
               by = c("stream", "replicate"))
  list(df = df_cor,
       r  = cor.test(df_cor$DLM_micro, df_cor$DLM_total, method = "pearson"))
}, error = function(e) NULL)

if (!is.null(cor_kz_dlm)) {
  r_val <- round(cor_kz_dlm$r$estimate, 3)
  p_val <- round(cor_kz_dlm$r$p.value,  3)
  message(glue("Kinzig DLM micro vs total: r = {r_val}, p = {p_val}"))

  p_dlm_cor <- ggplot(cor_kz_dlm$df, aes(DLM_micro, DLM_total)) +
    geom_point(alpha = 0.6, size = 2, colour = COL_RURAL) +
    geom_smooth(method = "lm", colour = COL_URBAN, se = TRUE, linewidth = 1) +
    annotate("text",
             x = min(cor_kz_dlm$df$DLM_micro) * 1.05,
             y = max(cor_kz_dlm$df$DLM_total) * 0.95,
             label = glue("r = {r_val},  p = {p_val}"),
             size = 3.5, hjust = 0) +
    labs(title = "DLM micro vs total correlation (Kinzig)",
         x = "DLM_micro (% °C⁻¹ d⁻¹)",
         y = "DLM_total (% °C⁻¹ d⁻¹)") +
    theme_pub()
  save_fig(p_dlm_cor, "Figure_DLM_micro_total_correlation.png")
}

## 29.4 Inter-replicate CV per stream ------------------------------------------
p_dlm_cv <- kinzig_micro_dlm %>%
  group_by(stream) %>%
  summarise(cv = sd(DLM, na.rm = TRUE) / mean(DLM, na.rm = TRUE) * 100,
            .groups = "drop") %>%
  mutate(stream = fct_reorder(stream, cv)) %>%
  ggplot(aes(stream, cv)) +
  geom_col(fill = COL_RURAL, alpha = 0.8) +
  coord_flip() +
  labs(title = "Inter-replicate CV per stream – Kinzig DLM_micro",
       x = NULL, y = "CV (%)") +
  theme_pub()
save_fig(p_dlm_cv, "Figure_DLM_CV_per_stream.png")


# =============================================================================
# 30  Variable preparation, transformation & collinearity  (Steps 3–5)
# =============================================================================
message("=== Section 30: Variable prep, transformation, collinearity ===")

## 30.1 Candidate predictors (no temperature) ----------------------------------
CAND_PRED <- c("pct_agri", "conductivity_u_s_cm",
               "no3_mg_l", "ortho_po4_mg_l",
               "cpom_pct", "fpom_pct", "psam_pct",
               "cd_mg_l", "cu_mg_l", "fe_mg_l",
               "ni_mg_l", "pb_mg_l", "zn_mg_l")

## 30.2 Helper: drop high-NA / high-zero predictors ---------------------------
filter_preds <- function(df, preds, na_thr = 0.40, zero_thr = 0.40) {
  preds <- intersect(preds, names(df))
  keep  <- sapply(preds, function(v) {
    x <- df[[v]]
    (mean(is.na(x)) < na_thr) & (mean(x == 0, na.rm = TRUE) < zero_thr)
  })
  preds[keep]
}

## 30.3 Helper: log10 transform if |skew| > 2 ----------------------------------
transform_df <- function(df, preds) {
  for (v in preds) {
    sk <- e1071::skewness(df[[v]], na.rm = TRUE)
    if (is.na(sk) || abs(sk) <= 2) next
    df[[v]] <- if (any(df[[v]] == 0, na.rm = TRUE)) log10(df[[v]] + 1) else log10(df[[v]])
    message(glue("  log10 transform: {v} (skew = {round(sk,2)})"))
  }
  # DLM response
  sk_dlm <- e1071::skewness(df$DLM, na.rm = TRUE)
  if (!is.na(sk_dlm) && abs(sk_dlm) > 2) {
    message(glue("  log10 transform: DLM (skew = {round(sk_dlm,2)})"))
    df$DLM <- log10(df$DLM)
  }
  df
}

## 30.4 Helper: pairwise |r|>0.7 then VIF>10 removal --------------------------
remove_collinear <- function(df, preds, r_thr = 0.7, vif_thr = 10) {
  preds <- intersect(preds, names(df))
  df_s  <- df %>% select(all_of(preds)) %>% drop_na()
  if (nrow(df_s) < 4 || ncol(df_s) < 2) return(preds)

  # Pairwise r
  cm <- cor(df_s, use = "complete.obs")
  high <- which(abs(cm) > r_thr & upper.tri(cm), arr.ind = TRUE)
  if (nrow(high) > 0) {
    dropped <- character(0)
    for (i in seq_len(nrow(high))) {
      v1 <- rownames(cm)[high[i,1]]; v2 <- colnames(cm)[high[i,2]]
      if (v1 %in% dropped || v2 %in% dropped) next
      mc1 <- mean(abs(cm[v1, setdiff(preds, v1)]), na.rm = TRUE)
      mc2 <- mean(abs(cm[v2, setdiff(preds, v2)]), na.rm = TRUE)
      drop_v <- if (mc1 >= mc2) v1 else v2
      message(glue("  |r|>{r_thr}: {v1} & {v2} → drop {drop_v}"))
      dropped <- c(dropped, drop_v)
    }
    preds <- setdiff(preds, dropped)
  }

  # VIF
  if (length(preds) >= 2) {
    df_vif <- df %>% select(DLM, all_of(preds)) %>% drop_na()
    if (nrow(df_vif) > length(preds) + 1) {
      repeat {
        m_v <- tryCatch(lm(as.formula(paste("DLM ~",
                            paste(preds, collapse = "+"))), data = df_vif),
                        error = function(e) NULL)
        if (is.null(m_v) || length(preds) < 2) break
        vf <- tryCatch(car::vif(m_v), error = function(e) NULL)
        if (is.null(vf) || max(vf, na.rm = TRUE) <= vif_thr) break
        drop_v <- names(which.max(vf))
        message(glue("  VIF>{vif_thr}: drop {drop_v} (VIF={round(max(vf),1)})"))
        preds <- setdiff(preds, drop_v)
      }
    }
  }
  preds
}

## 30.5 Apply to each DLM dataset ----------------------------------------------
# Kinzig micro
km_preds  <- filter_preds(kinzig_micro_dlm, CAND_PRED)
kinzig_micro_dlm_t <- transform_df(kinzig_micro_dlm, km_preds)
km_preds_f <- remove_collinear(kinzig_micro_dlm_t, km_preds)

# Kinzig total
kt_preds  <- filter_preds(kinzig_total_dlm, CAND_PRED)
kinzig_total_dlm_t <- transform_df(kinzig_total_dlm, kt_preds)
kt_preds_f <- remove_collinear(kinzig_total_dlm_t, kt_preds)

# Boye micro
bm_preds  <- filter_preds(boye_micro_dlm, CAND_PRED)
boye_micro_dlm_t <- transform_df(boye_micro_dlm, bm_preds)
bm_preds_f <- remove_collinear(boye_micro_dlm_t, bm_preds)

# Boye total
bt_preds  <- filter_preds(boye_total_dlm, CAND_PRED)
boye_total_dlm_t <- transform_df(boye_total_dlm, bt_preds)
bt_preds_f <- remove_collinear(boye_total_dlm_t, bt_preds)

message("Final predictors after collinearity filter:")
message("  Kinzig micro : ", paste(km_preds_f, collapse = ", "))
message("  Kinzig total : ", paste(kt_preds_f, collapse = ", "))
message("  Boye   micro : ", paste(bm_preds_f, collapse = ", "))
message("  Boye   total : ", paste(bt_preds_f, collapse = ", "))


# =============================================================================
# 31  Stability selection (Step 6)
# =============================================================================
message("=== Section 31: Stability selection (cutoff=0.70, PFER=1) ===")

run_stabsel_dlm <- function(df, preds, endpoint_nm) {
  df_s <- df %>% select(DLM, all_of(preds)) %>% drop_na()
  if (nrow(df_s) <= length(preds) + 2 || length(preds) < 2) {
    message(glue("  {endpoint_nm}: insufficient data, skip")); return(NULL)
  }
  X <- scale(as.matrix(df_s[, preds]))
  # drop zero-variance columns after scaling
  ok  <- apply(X, 2, function(v) sd(v, na.rm = TRUE) > 1e-10)
  X   <- X[, ok, drop = FALSE]
  if (ncol(X) < 2) { message(glue("  {endpoint_nm}: <2 valid preds")); return(NULL) }
  y   <- df_s$DLM

  stab <- tryCatch(
    stabs::stabsel(x = X, y = y,
                   fitfun    = stabs::glmnet.lasso,
                   cutoff    = 0.70,
                   PFER      = 1,
                   assumption = "unimodal"),
    error = function(e) { message("  stabsel error: ", e$message); NULL })
  if (is.null(stab)) return(NULL)

  stab_df <- tibble(
    predictor      = names(stab$max),
    selection_prob = as.numeric(stab$max),
    selected       = as.numeric(stab$max) >= 0.70,
    endpoint       = endpoint_nm
  ) %>% arrange(desc(selection_prob))

  p_stab <- ggplot(stab_df, aes(reorder(predictor, selection_prob),
                                 selection_prob, fill = selected)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0.70, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c(`FALSE` = "grey75", `TRUE` = COL_RURAL),
                      guide = "none") +
    coord_flip() + ylim(0, 1) +
    labs(title    = glue("Stability selection – DLM {endpoint_nm}"),
         subtitle = "Dashed = 0.70 threshold",
         x = NULL, y = "Selection probability") +
    theme_pub()

  save_fig(p_stab, glue("Figure_DLM_stabsel_{endpoint_nm}.png"), h = max(8, ncol(X) * 1.2))
  stab_df
}

stab_km <- run_stabsel_dlm(kinzig_micro_dlm_t, km_preds_f, "kinzig_micro")
stab_kt <- run_stabsel_dlm(kinzig_total_dlm_t, kt_preds_f, "kinzig_total")
stab_bm <- run_stabsel_dlm(boye_micro_dlm_t,   bm_preds_f, "boye_micro")
stab_bt <- run_stabsel_dlm(boye_total_dlm_t,   bt_preds_f, "boye_total")

stab_all <- bind_rows(stab_km, stab_kt, stab_bm, stab_bt)
if (nrow(stab_all) > 0) {
  write.csv(stab_all, "Table_DLM_stability_selection.csv", row.names = FALSE)
  message("Saved: Table_DLM_stability_selection.csv")
}


# =============================================================================
# 32  Full lm() models, marginal effects & DLM figures  (Steps 7–10)
# =============================================================================
message("=== Section 32: LM models, marginal effects & DLM figures ===")

## 32.1 Helper: get stability-selected predictors (fallback to top-3) ----------
sel_preds <- function(stab_df, fallback_preds) {
  if (is.null(stab_df)) return(head(fallback_preds, 3))
  sp <- stab_df %>% filter(selected) %>% pull(predictor)
  if (length(sp) == 0) return(head(fallback_preds, 3))
  sp
}

## 32.2 Helper: fit lm, export table, plot marginal effects -------------------
fit_and_plot_dlm <- function(df, preds, endpoint_nm) {
  sp   <- sel_preds(
    switch(endpoint_nm,
           kinzig_micro = stab_km, kinzig_total = stab_kt,
           boye_micro   = stab_bm, boye_total   = stab_bt, NULL),
    preds)
  sp   <- intersect(sp, names(df))
  if (length(sp) == 0) { message(glue("  {endpoint_nm}: no preds")); return(NULL) }

  df_m <- df %>%
    select(DLM, stream, all_of(sp)) %>% drop_na() %>%
    mutate(across(all_of(sp), ~scale(.)[, 1], .names = "{.col}_sc"))
  sc   <- paste0(sp, "_sc")

  if (nrow(df_m) <= length(sc) + 2) {
    message(glue("  {endpoint_nm}: too few rows")); return(NULL)
  }

  m <- lm(as.formula(paste("DLM ~", paste(sc, collapse = "+"))), data = df_m)
  message(glue("  {endpoint_nm}: R² = {round(summary(m)$r.squared, 3)}"))

  # Export coefficient table
  tidy_m <- broom::tidy(m, conf.int = TRUE) %>% mutate(endpoint = endpoint_nm)
  write.csv(tidy_m, glue("Table_DLM_lm_{endpoint_nm}.csv"), row.names = FALSE)

  # Marginal effect plots
  marg_plots <- lapply(sc, function(pr) {
    eff <- tryCatch(ggpredict(m, terms = paste0(pr, " [all]")),
                   error = function(e) NULL)
    if (is.null(eff)) return(NULL)
    ggplot(as.data.frame(eff), aes(x, predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = COL_RURAL, alpha = 0.15) +
      geom_line(colour = COL_RURAL, linewidth = 1.1) +
      labs(title = gsub("_sc$", "", pr),
           x = pr, y = "DLM (% °C⁻¹ d⁻¹)") +
      theme_pub(9)
  })
  marg_plots <- Filter(Negate(is.null), marg_plots)
  if (length(marg_plots) > 0) {
    p_marg <- wrap_plots(marg_plots, ncol = min(3, length(marg_plots))) +
      plot_annotation(title = glue("DLM marginal effects – {endpoint_nm}"))
    save_fig(p_marg, glue("Figure_DLM_marginal_{endpoint_nm}.png"),
             w = min(14 * length(marg_plots), 42), h = 12)
  }

  list(model = m, data = df_m, tidy = tidy_m, sc_preds = sc, name = endpoint_nm)
}

lm_km <- fit_and_plot_dlm(kinzig_micro_dlm_t, km_preds_f, "kinzig_micro")
lm_kt <- fit_and_plot_dlm(kinzig_total_dlm_t, kt_preds_f, "kinzig_total")
lm_bm <- fit_and_plot_dlm(boye_micro_dlm_t,   bm_preds_f, "boye_micro")
lm_bt <- fit_and_plot_dlm(boye_total_dlm_t,   bt_preds_f, "boye_total")

## 32.3 DLM similarity analysis (Step 9) --------------------------------------
get_betas <- function(lm_obj) {
  if (is.null(lm_obj)) return(tibble(term = character(), estimate = numeric()))
  broom::tidy(lm_obj$model) %>%
    filter(term != "(Intercept)") %>%
    mutate(term = gsub("_sc$", "", term)) %>%
    select(term, estimate)
}

beta_km_dlm  <- get_betas(lm_km) %>% rename(beta_micro = estimate)
beta_kt_dlm  <- get_betas(lm_kt) %>% rename(beta_total = estimate)
beta_kmac_dlm <- get_betas(
  tryCatch({
    df_mac <- kinzig_macro_dlm %>%
      select(DLM, stream, any_of(km_preds_f)) %>% drop_na() %>%
      mutate(across(any_of(km_preds_f), ~scale(.)[,1], .names = "{.col}_sc"))
    sc_mac <- paste0(intersect(km_preds_f, names(kinzig_macro_dlm)), "_sc")
    sc_mac <- intersect(sc_mac, names(df_mac))
    if (length(sc_mac) < 1 || nrow(df_mac) <= length(sc_mac) + 2) return(NULL)
    list(model = lm(as.formula(paste("DLM ~", paste(sc_mac, collapse="+"))),
                    data = df_mac))
  }, error = function(e) NULL)
) %>% rename(beta_macro = estimate)

dlm_sim <- beta_km_dlm %>%
  full_join(beta_kt_dlm,  by = "term") %>%
  full_join(beta_kmac_dlm, by = "term") %>%
  mutate(
    similarity_pct = apply(cbind(sign(beta_micro), sign(beta_total),
                                  sign(beta_macro)), 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      sum(combn(x, 2, FUN = function(p) p[1] == p[2])) /
        ncol(combn(x, 2)) * 100
    }),
    direction_class = case_when(
      similarity_pct == 100 ~ "Similar",
      similarity_pct ==   0 ~ "Opposing",
      TRUE                  ~ "Mixed")
  )

write.csv(dlm_sim, "Table_DLM_similarity.csv", row.names = FALSE)
message("Saved: Table_DLM_similarity.csv")

## 32.4 Combined coefficient comparison plot (Step 10) ------------------------
tidy_all_dlm <- bind_rows(
  if (!is.null(lm_km)) lm_km$tidy else NULL,
  if (!is.null(lm_kt)) lm_kt$tidy else NULL,
  if (!is.null(lm_bm)) lm_bm$tidy else NULL,
  if (!is.null(lm_bt)) lm_bt$tidy else NULL
) %>% filter(term != "(Intercept)") %>%
  mutate(term = gsub("_sc$", "", term))

if (nrow(tidy_all_dlm) > 0) {
  ep_cols <- c(kinzig_micro = COL_RURAL,   kinzig_total = "#66BB6A",
               boye_micro   = COL_URBAN,   boye_total   = "#EF9A9A")

  p_dlm_coef <- ggplot(tidy_all_dlm, aes(estimate, term, colour = endpoint)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.12, linewidth = 0.7,
                   position = position_dodge(0.4)) +
    geom_point(size = 3, position = position_dodge(0.4)) +
    scale_colour_manual(values = ep_cols, name = "Endpoint") +
    labs(title    = "DLM model coefficients across endpoints",
         subtitle = "Standardised predictors; error bars = 95 % CI",
         x = "Standardised coefficient", y = NULL) +
    theme_pub()
  save_fig(p_dlm_coef, "Figure_DLM_coef_comparison.png", w = 20, h = 14)
}

## 32.5 DLM similarity barplot ------------------------------------------------
if (nrow(dlm_sim) > 0 && any(!is.na(dlm_sim$similarity_pct))) {
  p_dlm_simbar <- dlm_sim %>%
    filter(!is.na(similarity_pct)) %>%
    mutate(term = fct_reorder(term, similarity_pct)) %>%
    ggplot(aes(term, similarity_pct, fill = direction_class)) +
    geom_col(width = 0.55) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c(Similar  = COL_RURAL,
                                  Opposing = COL_URBAN,
                                  Mixed    = COL_HIGH_FINE),
                      name = "Direction class") +
    coord_flip() + ylim(0, 100) +
    labs(title    = "DLM predictor similarity across endpoints",
         subtitle = "100 % = same direction in all DLM endpoints",
         x = NULL, y = "Similarity (%)") +
    theme_pub()
  save_fig(p_dlm_simbar, "Figure_DLM_similarity_barplot.png", h = 8)
}

## 32.6 DLM stream heatmap ----------------------------------------------------
dlm_stream_heat <- bind_rows(
  kinzig_micro_dlm %>% group_by(stream) %>%
    summarise(mean_DLM = mean(DLM, na.rm=TRUE), .groups="drop") %>%
    mutate(endpoint = "Kinzig DLM_micro"),
  kinzig_total_dlm %>% group_by(stream) %>%
    summarise(mean_DLM = mean(DLM, na.rm=TRUE), .groups="drop") %>%
    mutate(endpoint = "Kinzig DLM_total"),
  boye_micro_dlm %>% group_by(stream) %>%
    summarise(mean_DLM = mean(DLM, na.rm=TRUE), .groups="drop") %>%
    mutate(endpoint = "Boye DLM_micro"),
  boye_total_dlm %>% group_by(stream) %>%
    summarise(mean_DLM = mean(DLM, na.rm=TRUE), .groups="drop") %>%
    mutate(endpoint = "Boye DLM_total")
)

p_dlm_heat <- ggplot(dlm_stream_heat, aes(endpoint, stream, fill = mean_DLM)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "white", high = COL_RURAL, na.value = "grey88",
                      name = "Mean DLM\n(% °C⁻¹ d⁻¹)") +
  labs(title = "Mean DLM per stream and endpoint",
       x = "Endpoint", y = "Stream") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")
save_fig(p_dlm_heat, "Figure_DLM_stream_heatmap.png", w = 18, h = 20)

## 32.7 DLM violin by LandUse (Kinzig) and RestStatus (Boye) ------------------
p_dlm_lu <- ggplot(kinzig_micro_dlm,
                    aes(LandUse, DLM, fill = LandUse, colour = LandUse)) +
  geom_violin(alpha = 0.25, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA,
               colour = "grey30", fill = "white") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "white", colour = "grey20") +
  scale_fill_manual(values   = COL_PAL_LU, guide = "none") +
  scale_colour_manual(values = COL_PAL_LU, guide = "none") +
  labs(title = "DLM_micro by land use (Kinzig)",
       x = NULL, y = "DLM_micro (% °C⁻¹ d⁻¹)") +
  theme_pub()

p_dlm_rs <- ggplot(boye_micro_dlm,
                    aes(RestStatus, DLM, fill = RestStatus)) +
  geom_violin(alpha = 0.25, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA,
               colour = "grey30", fill = "white") +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "white", colour = "grey20") +
  scale_fill_manual(values = c(reference          = "#1976D2",
                                recently_restored  = "#F57C00",
                                long_restored      = "#388E3C"),
                    guide = "none") +
  labs(title = "DLM_micro by restoration status (Boye)",
       x = NULL, y = "DLM_micro (% °C⁻¹ d⁻¹)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

panel_dlm_violin <- (p_dlm_lu | p_dlm_rs) +
  plot_annotation(tag_levels = "A",
                  title = "DLM_micro by catchment grouping")
save_fig(panel_dlm_violin, "Figure_DLM_violin_LU_RS.png", w = 28, h = 12)


# =============================================================================
# 33  Summary of all output files
# =============================================================================
message("\n=== Analysis complete ===")
all_out <- list.files(pattern = "^(Figure|Table|Stabsel).*\\.(png|csv)$")
message(glue("Total output files in working directory: {length(all_out)}"))
cat(paste(" •", sort(all_out)), sep = "\n")
