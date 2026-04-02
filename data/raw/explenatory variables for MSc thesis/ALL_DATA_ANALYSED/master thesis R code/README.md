# Master Thesis R Code – Agube Ikenna Stephen
**Project:** RESIST Field Study – Leaf Decomposition in Urban Streams (Kinzig & Boye, 2025)
**Institution:** University of Duisburg-Essen

---

## Folder Structure

### R_scripts/
All R and RMarkdown scripts for the thesis analysis, numbered in order of use.

| File | What it does |
|---|---|
| `00_session_info_packages_used.txt` | Lists all R packages and versions used – essential for reproducibility |
| `01_kinzig_decomposition_analysis_FINAL.R` | The final, clean analysis script for leaf decomposition rates at Kinzig sites |
| `02_kinzig_decomposition_analysis_v1.R` | First version of the Kinzig analysis – useful for comparison |
| `03_kinzig_decomposition_analysis_v2.R` | Second version with updated model structure |
| `04_transport_controls_check.R` | Script to check and process transport control leaf bags |
| `05_standard_vs_natural_comparison.Rmd` | RMarkdown comparing standard (cotton strip / DecotaB) vs natural (leaf) decomposition |
| `06_general_analysis.R` | General exploratory analysis script – combines Boye and Kinzig for data prep |
| `07_boye_decomposition_analysis.R` | **Boye-specific decomposition analysis** – uses RESTORATION STATUS as predictor (reference vs recently restored vs long-restored), NOT land use |

---

### data/master_data/
The raw and combined field master datasets. Start here when building any analysis.

| File | What it contains |
|---|---|
| `RESIST_field_master_data_BoyeKinzig_all_chemistry.xlsx` | Full master data combining Boye and Kinzig sites with all chemical variables |
| `RESIST_kinzig_chemistry_data_2025.xlsx` | Physicochemical data for Kinzig sites specifically |
| `RESIST_field_raw_data_2025.xlsx` | Unprocessed raw field data from the 2025 sampling campaign |

---

### data/processed/
Cleaned and processed data files ready for analysis in R.

| File | What it contains |
|---|---|
| `field_initial_leaf_mass_2025.xlsx` | Starting dry weights of leaf bags before deployment |
| `field_processed_samples_2025.xlsx` | Processed sample data from leaf bag retrieval |
| `field_transport_controls_2025.xlsx` | Transport control bag data (accounting for handling/drying loss) |
| `field_invertebrates_coarse_leaf_bags_2025.xlsx` | Macroinvertebrate counts from coarse-mesh leaf bags |
| `field_site_information_kinzig_2025.xlsx` | Site coordinates, catchment details, and reach characteristics |
| `physicochemical_data_kinzig_2025.xlsx` | Water chemistry measurements (pH, conductivity, oxygen, etc.) |
| `substrate_estimation_kinzig_2025_clean.xlsx` | Cleaned substrate composition estimates at Kinzig sites |
| `substrate_estimation_boye_clean.xlsx` | Cleaned substrate composition estimates at Boye sites |
| `macroinvertebrate_boye_2025_processed.xlsx` | Processed macroinvertebrate data from Boye sites |
| `RESIST_field_variables_master_march2025.xlsx` | Master variable reference with all environmental predictors |
| `field_date_summary_boye_kinzig.xlsx` | Summary of field deployment and retrieval dates |
| `kinzig_processed_cleaned.csv` | Final cleaned dataset used in the mixed models (CSV format for R) |
| `kinzig_decomposition_rates_summary.csv` | Calculated decomposition rates (k values) per site and leaf type |
| `boye_processed_samples_mar2025.xlsx` | Boye water chemistry data (pH, conductivity, DO, metals, nutrients) |
| `boye_temperature_logger_march2025.csv` | Temperature logger data for Boye – use this to calculate degree-days accurately |
| `water_samples_analysis_march2025.xlsx` | Water sample analysis results for March 2025 (Boye + Kinzig) |
| `land_use_summary_2022_BoyeKinzig.xlsx` | Land use data for both catchments (background variable) |
| `SFB_combined_data_BoyeKinzig.xlsx` | Combined dataset across both catchments for cross-catchment comparisons |
| `kinzig_processed_samples_mar2025.xlsx` | Kinzig processed samples final file for March 2025 |
| `site_basic_data_BoyeKinzig_2025.xlsx` | Sampling site coordinates, restoration year, and position for all Boye and Kinzig sites |

---

### thesis_documents/
Core writing and planning documents.

| File | What it is |
|---|---|
| `thesis_registration_form.pdf` | Official MSc thesis registration form |
| `university_master_thesis_guidelines.pdf` | UDE formatting and submission guidelines – read before writing |
| `thesis_proposal_draft.docx` | Initial thesis proposal |
| `thesis_introduction_and_methods.docx` | Draft introduction and methods chapters |
| `thesis_hypothesis_original.docx` | First version of the thesis hypotheses |
| `thesis_hypothesis_revised.docx` | Revised, updated hypotheses |
| `kinzig_catchment_hypothesis_final.docx` | Final hypotheses specific to the Kinzig catchment study |
| `variable_reference_guide.docx` | Explains what each explanatory variable is and how it was measured |
| `leaf_decomposition_methods_guide.docx` | Step-by-step methods guide for leaf decomposition processing |
| `data_cleaning_rules.docx` | Rules applied during data cleaning – document any deviations from this |

---

### literature/key_articles/
The core papers for understanding and citing your thesis work.

**Leaf decomposition methods and theory**
- `Gessner_2002_litter_breakdown_functional_stream_integrity.pdf` – The foundational argument for using leaf breakdown rates as a functional ecosystem measure. Cite in your introduction.
- `Gessner_1997_fungal_biomass_sporulation_leaf_decomposition.pdf` – Explains fungal contribution to decomposition, relevant to your microbial decomposition rate.
- `Barlocher_1985_fungi_stream_invertebrate_nutrition.pdf` – Classic paper on the role of fungi in making leaf litter nutritious for invertebrates.
- `Hunting_2015_DecotaB_standard_substrate_microbial_decomposition.pdf` – Describes the DecotaB method used in the Standard vs Natural comparison.

**Standard vs Natural decomposition substrates**
- `Schreiner_2023_Standard_vs_Natural_leaf_decomposition.pdf` – Directly related to your study design. Compares cotton strips and DecotaB against natural leaves.
- `Schreiner_2018_microbial_recovery_urbanized_streams.pdf` – Shows how long microbial functions take to recover in restored urban streams.

**Multiple stressors and urban stream ecology**
- `Medina_Madariaga_2024_multiple_stressor_leaf_litter_decomposition.pdf` – Recent paper on multiple stressor effects on decomposition. Essential for your discussion.
- `Multiple_stressors_microbial_decomposer_litter_breakdown.pdf` – Covers stressor interactions affecting decomposer communities.
- `Madge_Pimentel_2024_urban_stream_ecosystem_response.pdf` – Urban stream ecosystem functioning under stress.
- `Urban_Stream_Syndrome_impact_on_macroinvertebrates.pdf` – Evaluates how urban stream degradation affects macroinvertebrate communities.

**Land use, catchment, and water quality**
- `Land_use_catchment_effect_water_quality_plant_species.pdf` – How catchment land use shapes water quality and biotic communities.
- `Stream_macroinvertebrates_restored_vs_impacted_catchments_climate_landuse.pdf` – Shows how macroinvertebrate communities respond to climate, land use, and flow across restored vs degraded sites.
- `Link_2017_dilution_factors_german_wastewater_treatment.pdf` – German-specific paper on WWTP dilution factors. Relevant to your Kinzig sites downstream of WWTPs.

**Parasites and stressor cascades**
- `Stream_degradation_recovery_parasite_communities_multiple_stressors.pdf` – Useful for thinking about indirect ecological effects of stream degradation.
- `Fernandez_2015_fungicides_decomposer_communities_leaf_breakdown.pdf` – Fungicide effects on decomposer communities – relevant if your sites receive agricultural runoff.

**Biogeography and global context**
- `Zhang_2019_Global_Ecology_Biogeography.pdf` – Broad-scale context for stream ecology patterns.
- `Environmental_drivers_leaf_breakdown_rate_urban_watershed.pdf` – Identifies which environmental variables best predict leaf breakdown in urban systems.

---

### literature/msc_theses/
MSc theses to read for structure, framing, and methods inspiration.

| File | Why it is useful |
|---|---|
| `Stefansdottir_2010_MScThesis_IMPORTANT.pdf` | Flagged as important – read this one first |
| `Icaza_MScThesis_leaf_decomposition_urban_streams.pdf` | Directly relevant: leaf decomposition in urban streams |
| `Brokaw_2025_MScThesis_stream_ecology.pdf` | Recent stream ecology thesis – useful for discussion framing |
| `Hashmi_MScThesis_freshwater_ecology.pdf` | Freshwater ecology methods and writing structure |
| `Stragnefeldt_MScThesis_stream_ecology.pdf` | Stream ecology thesis – study design and methods reference |
| `VT2025_ES2530_NK_MScThesis.pdf` | Recent environmental science MSc thesis |

---

## How to start your R analysis

1. Read `00_session_info_packages_used.txt` and install any packages you are missing.
2. Open `01_kinzig_decomposition_analysis_FINAL.R` – this is the most current working script.
3. The script reads from `data/processed/kinzig_processed_cleaned.csv` and `data/master_data/`.
4. Cross-reference `thesis_documents/variable_reference_guide.docx` for any variable you do not recognise.
5. If you need to recheck transport controls, run `04_transport_controls_check.R` first.

---

*Last updated: March 2026*
