# ================================================================================================
# CLIF ARF National Estimation | Federated Site Export Script (End-to-End)
# PI: Peter Graffy
#
# Outputs (PHI-safe):
#   1) site_county_year.csv            : site_name × county_fips × year counts + missingness + coverage
#   2) site_county_year_age_sex.csv    : site_name × county_fips × year × age_band × sex counts
#   3) site_qc_summary.csv             : QC tallies for phenotype + missingness + data availability
#
# Requires CLIF v2.1+ tables:
#   clif_patient, clif_hospitalization, clif_adt, clif_vitals, clif_labs, clif_respiratory_support
#
# Notes:
# - County must be county-of-residence (preferred). Script uses hospitalization.county_code by default.
# - ARF phenotype is physiologic (SpO2/FIO2, PaO2/FIO2, PaCO2/pH) within ±24h of first ICU in.
# - Mortality is in-hospital death time using patient.death_dttm with discharge fallback to last vital.
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(readr)
  library(arrow)
  library(fst)
  library(data.table)
  library(glue)
  library(janitor)
})

# ------------------------------- 0) Config -------------------------------------------------------

# Load configuration utility
source("utils/config.R")
repo <- config$repo
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type

print(paste("Site Name:", site_name))
print(paste("Tables Path:", tables_path))
print(paste("File Type:", file_type))

# Study window (as in SAP; adjust end year to 2025 if needed)
START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE   <- as.POSIXct("2025-12-31 23:59:59", tz = "UTC")

# ARF phenotype parameters (from your prior script)
WINDOW_H        <- 24     # ± hours around ICU in
N_SPO2_MIN      <- 6      # continuous SpO2 proxy
ROOM_AIR_FIO2   <- 0.21
JOIN_NEAR_H     <- 1      # SpO2/PaO2 ↔ FiO2 pairing window in hours
HYPERPAIR_H     <- 2      # PaCO2 ↔ pH pairing window in hours
MIN_ICU_LOS_H   <- 24     # ICU LOS threshold for inclusion

# If TRUE, enforce demo presence (age, sex, race) and geo presence before counting A/Y/D
ENFORCE_DEMO_GEO_FILTERS <- TRUE

out_dir <- file.path("output", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
export_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_site_county_year      <- file.path(out_dir, glue("site_county_year_{site_name}_{export_stamp}.csv"))
out_site_county_year_as   <- file.path(out_dir, glue("site_county_year_age_sex_{site_name}_{export_stamp}.csv"))
out_site_qc               <- file.path(out_dir, glue("site_qc_summary_{site_name}_{export_stamp}.csv"))

# ------------------------------- 1) IO helpers ---------------------------------------------------

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext, " for path: ", path))
}

# Robust datetime parser (fast; avoids locale issues)
safe_ts <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = tz))
  if (is.numeric(x)) {
    x2 <- ifelse(x > 1e12, x/1000, x)  # ms vs sec
    return(as.POSIXct(x2, origin = "1970-01-01", tz = tz))
  }
  # try fast parse
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c("ymd_HMS","ymd_HM","ymd","ymdTz","ymdT","mdy_HMS","mdy_HM","mdy","dmy_HMS","dmy_HM","dmy","HMS"),
    tz = tz, quiet = TRUE
  ))
}

safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

get_tbl_from_dir <- function(tbl_base) {
  # Accept clif_patient or patient; normalize to clif_* naming
  wanted <- tolower(tbl_base)
  if (!startsWith(wanted, "clif_")) wanted <- paste0("clif_", wanted)
  
  files <- list.files(tables_path, full.names = TRUE, recursive = TRUE)
  files <- files[grepl("\\.(csv|parquet|fst)$", files, ignore.case = TRUE)]
  
  bn <- tolower(tools::file_path_sans_ext(basename(files)))
  # normalize "patient" -> "clif_patient"
  bn_norm <- ifelse(startsWith(bn, "clif_"), bn, paste0("clif_", bn))
  
  hit <- files[bn_norm == wanted]
  if (length(hit) != 1) {
    stop(glue("Could not uniquely locate {wanted} in {tables_path}. Matches: {length(hit)}"))
  }
  janitor::clean_names(read_any(hit))
}

normalize_county_fips <- function(x) {
  # Expect 5-digit numeric string; anything else -> NA
  x <- as.character(x)
  x <- str_replace_all(x, "[^0-9]", "")
  x <- ifelse(nchar(x) == 5, x, NA_character_)
  x
}

age_band_4 <- function(age) {
  cut(age,
      breaks = c(18, 40, 65, 75, Inf),
      right = FALSE,
      labels = c("18-39", "40-64", "65-74", "75+"))
}

# ------------------------------- 2) Load minimal CLIF tables ------------------------------------

patient <- get_tbl_from_dir("clif_patient") %>%
  transmute(
    patient_id,
    birth_date = safe_date(birth_date),
    sex_category,
    race_category,
    ethnicity_category,
    death_dttm
  ) %>%
  mutate(death_dttm = safe_ts(death_dttm))

hospitalization <- get_tbl_from_dir("clif_hospitalization") %>%
  transmute(
    patient_id,
    hospitalization_id,
    admission_dttm  = safe_ts(admission_dttm),
    discharge_dttm  = safe_ts(discharge_dttm),
    age_at_admission,
    discharge_category,
    county_code,                # expected to reflect residence county; site-specific if not
    census_tract,
    zipcode_nine_digit,
    zipcode_five_digit
  )

adt <- get_tbl_from_dir("clif_adt") %>%
  transmute(
    hospitalization_id,
    in_dttm  = safe_ts(in_dttm),
    out_dttm = safe_ts(out_dttm),
    location_category
  )

vitals <- get_tbl_from_dir("clif_vitals") %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    vital_category = tolower(vital_category),
    vital_value
  ) %>%
  mutate(vital_value = suppressWarnings(as.numeric(vital_value)))

labs <- get_tbl_from_dir("clif_labs") %>%
  transmute(
    hospitalization_id,
    lab_result_dttm = safe_ts(lab_result_dttm),
    lab_category = tolower(lab_category),
    lab_value_numeric
  ) %>%
  mutate(lab_value_numeric = suppressWarnings(as.numeric(lab_value_numeric)))

resp_support <- get_tbl_from_dir("clif_respiratory_support") %>%
  transmute(
    hospitalization_id,
    recorded_dttm = safe_ts(recorded_dttm),
    fio2_set
  ) %>%
  mutate(fio2_set = suppressWarnings(as.numeric(fio2_set)))

# ------------------------------- 3) ICU index time per hospitalization ----------------------------

icu_segments <- adt %>%
  mutate(loccat = tolower(coalesce(location_category, "")),
         is_icu = str_detect(loccat, "icu")) %>%
  filter(is_icu) %>%
  filter(!is.na(in_dttm), !is.na(out_dttm), out_dttm > in_dttm)

icu_bounds <- icu_segments %>%
  group_by(hospitalization_id) %>%
  summarise(
    first_icu_in = suppressWarnings(min(in_dttm, na.rm = TRUE)),
    last_icu_out = suppressWarnings(max(out_dttm, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_icu_in = as.POSIXct(first_icu_in, tz="UTC"),
    last_icu_out = as.POSIXct(last_icu_out, tz="UTC"),
    icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units = "hours"))
  )

# ------------------------------- 4) Base candidate hospitalizations -------------------------------

base <- hospitalization %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
  left_join(patient %>% dplyr::select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
            by = "patient_id") %>%
  mutate(
    age_years = coalesce(
      suppressWarnings(as.numeric(age_at_admission)),
      ifelse(!is.na(birth_date),
             as.numeric(floor((as.Date(admission_dttm) - birth_date)/365.25)), NA_real_)
    ),
    adult = !is.na(age_years) & age_years >= 18,
    has_demo = !(is.na(age_years) | is.na(sex_category) | is.na(race_category)),
    has_geo = !is.na(census_tract) | !is.na(zipcode_nine_digit) | !is.na(zipcode_five_digit) | !is.na(county_code),
    icu_24h = !is.na(icu_los_hours) & icu_los_hours >= MIN_ICU_LOS_H
  ) %>%
  filter(!is.na(first_icu_in),
         first_icu_in >= START_DATE,
         first_icu_in <= END_DATE)

# Window boundaries
win <- base %>%
  transmute(
    hospitalization_id,
    win_start = first_icu_in - dhours(WINDOW_H),
    win_end   = first_icu_in + dhours(WINDOW_H)
  )

# ------------------------------- 5) Windowed signals --------------------------------------------

vitals_win <- vitals %>%
  filter(vital_category == "spo2") %>%
  inner_join(win, by = "hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  dplyr::select(hospitalization_id, recorded_dttm, spo2 = vital_value)

labs_win <- labs %>%
  filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
  inner_join(win, by = "hospitalization_id") %>%
  filter(lab_result_dttm >= win_start, lab_result_dttm <= win_end) %>%
  transmute(hospitalization_id, lab_result_dttm, lab_category, val = lab_value_numeric)

fio2_win <- resp_support %>%
  inner_join(win, by = "hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  dplyr::select(hospitalization_id, recorded_dttm, fio2_set)

# ------------------------------- 6) Pairing logic (data.table for speed) --------------------------

setDT(vitals_win); setDT(fio2_win); setDT(labs_win)

# Ensure POSIXct for rolling joins
vitals_win[, recorded_dttm := as.POSIXct(recorded_dttm, tz="UTC")]
fio2_win[, recorded_dttm   := as.POSIXct(recorded_dttm, tz="UTC")]
labs_win[, lab_result_dttm := as.POSIXct(lab_result_dttm, tz="UTC")]

# SpO2 ↔ FiO2 nearest within JOIN_NEAR_H
vdt <- vitals_win[, .(hospitalization_id, spo2_time = recorded_dttm, spo2 = spo2)]
fdt <- fio2_win[, .(hospitalization_id, fio2_time = recorded_dttm, fio2_set = fio2_set)]

setkey(vdt, hospitalization_id, spo2_time)
setkey(fdt, hospitalization_id, fio2_time)

vdt[, spo2_time_keep := spo2_time]

spo2_fio2 <- fdt[
  vdt, roll = "nearest", on = .(hospitalization_id, fio2_time = spo2_time), nomatch = 0L
][
  , timediff_h := abs(as.numeric(difftime(spo2_time_keep, fio2_time, units="hours")))
][
  timediff_h <= as.numeric(JOIN_NEAR_H),
  .(hospitalization_id, spo2_time = spo2_time_keep, fio2_time, spo2, fio2_set,
    on_room_air = !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6,
    timediff_h)
]

# PaO2 ↔ FiO2 nearest within JOIN_NEAR_H
po2dt <- labs_win[lab_category == "po2_arterial",
                  .(hospitalization_id, po2_time = lab_result_dttm, po2 = val)]
setkey(po2dt, hospitalization_id, po2_time)
po2dt[, po2_time_keep := po2_time]

po2_fio2 <- fdt[
  po2dt, roll = "nearest", on = .(hospitalization_id, fio2_time = po2_time), nomatch = 0L
][
  , timediff_h := abs(as.numeric(difftime(po2_time_keep, fio2_time, units="hours")))
][
  timediff_h <= as.numeric(JOIN_NEAR_H),
  .(hospitalization_id, po2_time = po2_time_keep, fio2_time, po2, fio2_set,
    pf_ratio = fifelse(!is.na(fio2_set) & fio2_set > 0, po2 / fio2_set, as.numeric(NA)),
    timediff_h)
]

# PaCO2 ↔ pH nearest within HYPERPAIR_H
pco2dt <- labs_win[lab_category == "pco2_arterial",
                   .(hospitalization_id, pco2_time = lab_result_dttm, pco2 = val)]
phdt   <- labs_win[lab_category == "ph_arterial",
                   .(hospitalization_id, ph_time = lab_result_dttm, ph = val)]

setkey(pco2dt, hospitalization_id, pco2_time)
setkey(phdt,   hospitalization_id, ph_time)

pco2dt[, pco2_time_keep := pco2_time]
phdt[,   ph_time_keep   := ph_time]

hyper_pairs <- phdt[
  pco2dt, roll="nearest", on=.(hospitalization_id, ph_time = pco2_time), nomatch = 0L
][
  , timediff_h := abs(as.numeric(difftime(pco2_time_keep, ph_time, units="hours")))
][
  timediff_h <= as.numeric(HYPERPAIR_H),
  .(hospitalization_id,
    pco2_time = pco2_time_keep, ph_time = ph_time_keep,
    pco2, ph,
    hyper_pair_hit = (pco2 >= 45 & ph < 7.35),
    timediff_h)
]

# Convert back to tibbles for dplyr summaries
spo2_fio2 <- as_tibble(spo2_fio2)
po2_fio2  <- as_tibble(po2_fio2)
hyper_pairs <- as_tibble(hyper_pairs)

# ------------------------------- 7) ARF flags per hospitalization --------------------------------

hypox_roomair_spo2 <- spo2_fio2 %>%
  mutate(hit = (spo2 < 90 & on_room_air)) %>%
  group_by(hospitalization_id) %>%
  summarise(any_spo2_roomair_hit = any(hit, na.rm = TRUE),
            spo2_n_pairs = n(),
            .groups = "drop")

hypox_roomair_po2 <- po2_fio2 %>%
  mutate(hit = (po2 <= 60 & !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6)) %>%
  group_by(hospitalization_id) %>%
  summarise(any_po2_roomair_hit = any(hit, na.rm = TRUE),
            .groups = "drop")

hypox_pf <- po2_fio2 %>%
  mutate(hit = (!is.na(pf_ratio) & pf_ratio <= 300)) %>%
  group_by(hospitalization_id) %>%
  summarise(any_pf_hit = any(hit, na.rm = TRUE), .groups = "drop")

hyper_flags <- hyper_pairs %>%
  group_by(hospitalization_id) %>%
  summarise(any_hyper_pair = any(hyper_pair_hit, na.rm = TRUE), .groups = "drop")

# Data availability: ABG OR continuous SpO2
abg_avail <- labs_win %>%
  filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
  distinct(hospitalization_id) %>%
  mutate(has_abg = TRUE)

spo2_density <- as_tibble(vitals_win) %>%
  group_by(hospitalization_id) %>%
  summarise(n_spo2 = n(), .groups = "drop") %>%
  mutate(has_cont_spo2 = n_spo2 >= N_SPO2_MIN)

data_avail <- base %>%
  dplyr::select(hospitalization_id) %>%
  left_join(abg_avail, by = "hospitalization_id") %>%
  left_join(spo2_density, by = "hospitalization_id") %>%
  mutate(
    has_abg = coalesce(has_abg, FALSE),
    has_cont_spo2 = coalesce(has_cont_spo2, FALSE),
    meets_data_rule = has_abg | has_cont_spo2
  )

flags <- base %>%
  left_join(hypox_roomair_spo2, by = "hospitalization_id") %>%
  left_join(hypox_roomair_po2,  by = "hospitalization_id") %>%
  left_join(hypox_pf,           by = "hospitalization_id") %>%
  left_join(hyper_flags,        by = "hospitalization_id") %>%
  left_join(data_avail %>% dplyr::select(hospitalization_id, meets_data_rule, has_abg, has_cont_spo2), by = "hospitalization_id") %>%
  mutate(
    any_hypox = coalesce(any_spo2_roomair_hit, FALSE) |
      coalesce(any_po2_roomair_hit, FALSE) |
      coalesce(any_pf_hit, FALSE),
    any_hypercap = coalesce(any_hyper_pair, FALSE),
    arf_criterion_met = any_hypox | any_hypercap,
    mixed_arf = any_hypox & any_hypercap,
    year = year(admission_dttm),
    age_band = age_band_4(age_years),
    county_fips = normalize_county_fips(county_code),
    missing_county = is.na(county_fips)
  )

# Inclusion mask (this defines eligible ICU denominator A_ct)
# --- County FIPS is required for denominators/outcomes (county-of-residence analysis) ---
flags <- flags %>%
  mutate(
    county_fips = normalize_county_fips(county_code),
    has_county  = !is.na(county_fips),
    
    # ICU base inclusion for denominator: adult + ICU LOS>=24h + county present
    include_all_icu = adult & icu_24h & has_county,
    
    # Phenotype eligibility: data rule only (ABG OR continuous SpO2) within ±24h
    include_elig = include_all_icu & meets_data_rule,
    
    # ARF case: physiologic criteria met
    include_arf  = include_elig & arf_criterion_met
  )

# ------------------------------- 8) Mortality definition (your logic) ----------------------------

# last vital per hospitalization for discharge fallback
vitals_last <- vitals %>%
  filter(hospitalization_id %in% flags$hospitalization_id) %>%
  mutate(vital_recorded_ts = safe_ts(recorded_dttm)) %>%
  filter(!is.na(vital_recorded_ts)) %>%
  group_by(hospitalization_id) %>%
  summarise(last_vital_dttm = max(vital_recorded_ts), .groups = "drop")

# Build final death timestamp per hospitalization
final_outcome_times <- hospitalization %>%
  filter(hospitalization_id %in% flags$hospitalization_id) %>%
  transmute(
    patient_id,
    hospitalization_id,
    discharge_category,
    discharge_dttm = safe_ts(discharge_dttm),
    discharge_cat_low = tolower(discharge_category)
  ) %>%
  left_join(patient %>% dplyr::select(patient_id, death_dttm), by = "patient_id") %>%
  left_join(vitals_last, by = "hospitalization_id") %>%
  mutate(
    death_dttm_final = case_when(
      discharge_cat_low %in% c("expired","hospice") & is.na(death_dttm) ~ last_vital_dttm,
      TRUE ~ death_dttm
    ),
    death_source = case_when(
      !is.na(death_dttm) ~ "patient_death_dttm",
      discharge_cat_low %in% c("expired","hospice") & !is.na(last_vital_dttm) ~ "discharge_fallback_last_vital",
      discharge_cat_low %in% c("expired","hospice") &  is.na(last_vital_dttm) ~ "discharge_fallback_missing_last_vital",
      TRUE ~ "no_death"
    ),
    death_ts = safe_ts(death_dttm_final)
  )

# Join mortality to flags; compute in-hospital and 30d
flags <- flags %>%
  left_join(final_outcome_times %>% dplyr::select(hospitalization_id, death_ts, death_source), by = "hospitalization_id") %>%
  mutate(
    in_hosp_death = as.integer(!is.na(death_ts) & death_ts >= admission_dttm & death_ts <= discharge_dttm),
    death_30d     = as.integer(!is.na(death_ts) & death_ts <= (admission_dttm + days(30)))
  )

# ------------------------------- 9) Aggregations for export --------------------------------------

# Coverage metadata (year-level)
coverage_year <- flags %>%
  filter(include_all_icu) %>%
  mutate(date = as.Date(admission_dttm)) %>%
  group_by(year) %>%
  summarise(
    first_date = min(date, na.rm = TRUE),
    last_date  = max(date, na.rm = TRUE),
    coverage_days = as.integer(difftime(last_date, first_date, units="days")) + 1L,
    .groups = "drop"
  )

# Primary county×year table
site_county_year <- flags %>%
  group_by(county_fips, year) %>%
  summarise(
    A_all_ct  = sum(include_all_icu, na.rm = TRUE),
    A_elig_ct = sum(include_elig, na.rm = TRUE),
    Y_ct      = sum(include_arf, na.rm = TRUE),
    D_ct      = sum(include_arf & in_hosp_death == 1L, na.rm = TRUE),
    D30_ct    = sum(include_arf & death_30d == 1L, na.rm = TRUE),
    
    # QC-only (should be zero if has_county enforced in include_* flags)
    A_all_missing_county  = sum(include_all_icu & !has_county, na.rm = TRUE),
    A_elig_missing_county = sum(include_elig & !has_county, na.rm = TRUE),
    Y_missing_county      = sum(include_arf & !has_county, na.rm = TRUE),
    D_missing_county      = sum(include_arf & in_hosp_death == 1L & !has_county, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(site_name = site_name) %>%
  left_join(coverage_year, by = "year") %>%
  arrange(year, county_fips)

# Recommended stratified county×year×age×sex
site_county_year_age_sex <- flags %>%
  mutate(sex = as.character(sex_category),
         age_band = as.character(age_band)) %>%
  group_by(county_fips, year, age_band, sex) %>%
  summarise(
    A_all_ctg  = sum(include_all_icu, na.rm = TRUE),
    A_elig_ctg = sum(include_elig, na.rm = TRUE),
    Y_ctg      = sum(include_arf, na.rm = TRUE),
    D_ctg      = sum(include_arf & in_hosp_death == 1L, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site_name = site_name) %>%
  arrange(year, county_fips, age_band, sex)

# ------------------------------- 10) QC summary --------------------------------------------------

qc <- tibble(
  site_id = site_name,
  start_date = as.character(START_DATE),
  end_date   = as.character(END_DATE),
  
  n_candidates = nrow(base),
  
  # Denominators + outcomes
  n_all_icu_geo = sum(flags$include_all_icu, na.rm = TRUE),
  n_elig_icu    = sum(flags$include_elig, na.rm = TRUE),
  n_arf         = sum(flags$include_arf, na.rm = TRUE),
  
  # Useful ratios
  prop_elig_among_all = ifelse(sum(flags$include_all_icu, na.rm = TRUE) > 0,
                               sum(flags$include_elig, na.rm = TRUE) / sum(flags$include_all_icu, na.rm = TRUE),
                               NA_real_),
  prop_arf_among_elig = ifelse(sum(flags$include_elig, na.rm = TRUE) > 0,
                               sum(flags$include_arf, na.rm = TRUE) / sum(flags$include_elig, na.rm = TRUE),
                               NA_real_),
  
  # ARF subtype counts among ARF
  n_hypoxemic   = sum(flags$include_arf & flags$any_hypox, na.rm = TRUE),
  n_hypercapnic = sum(flags$include_arf & flags$any_hypercap, na.rm = TRUE),
  n_mixed       = sum(flags$include_arf & flags$mixed_arf, na.rm = TRUE),
  
  # Data availability mix (computed among ALL ICU with geo)
  n_has_abg = sum(flags$include_all_icu & flags$has_abg, na.rm = TRUE),
  n_has_cont_spo2 = sum(flags$include_all_icu & flags$has_cont_spo2, na.rm = TRUE),
  n_has_abg_or_cont_spo2 = sum(flags$include_all_icu & flags$meets_data_rule, na.rm = TRUE),
  
  # Pairing volumes (helps detect FiO2 sparsity)
  n_spo2_fio2_pairs = nrow(spo2_fio2),
  n_po2_fio2_pairs  = nrow(po2_fio2),
  n_hyper_pairs     = nrow(hyper_pairs),
  
  # Geography missingness (QC only; should be ~0 because include_all_icu requires county)
  n_missing_county_all_icu = sum(flags$include_all_icu & flags$missing_county, na.rm = TRUE),
  
  # Mortality provenance among ARF
  n_arf_deaths_in_hosp = sum(flags$include_arf & flags$in_hosp_death == 1L, na.rm = TRUE),
  n_death_source_patient = sum(flags$include_arf & flags$death_source == "patient_death_dttm", na.rm = TRUE),
  n_death_source_fallback = sum(flags$include_arf & flags$death_source == "discharge_fallback_last_vital", na.rm = TRUE),
  n_death_fallback_missing_last_vital = sum(flags$include_arf & flags$death_source == "discharge_fallback_missing_last_vital", na.rm = TRUE)
)

qc <- qc %>%
  mutate(
    n_all_icu_geo = sum(flags$include_all_icu, na.rm = TRUE),
    n_elig_icu    = sum(flags$include_elig, na.rm = TRUE),
    n_arf         = sum(flags$include_arf, na.rm = TRUE)
  )

# ------------------------------- 11) Write outputs -----------------------------------------------

# Important: these exports are aggregated only; no IDs.
readr::write_csv(site_county_year, out_site_county_year, na = "")
readr::write_csv(site_county_year_age_sex, out_site_county_year_as, na = "")
readr::write_csv(qc, out_site_qc, na = "")

message("✅ Done. Wrote PHI-safe exports:")
message("  - ", out_site_county_year)
message("  - ", out_site_county_year_as)
message("  - ", out_site_qc)

# ------------------------------- 12) Optional: quick console sanity ------------------------------

message(glue("Site {site_name}: base included ICU hosp = {qc$n_elig_icu}, ARF = {qc$n_arf}, in-hosp ARF deaths = {qc$n_arf_deaths_in_hosp}"))

# =========================
# CONSORT-style flow export
# =========================

# Helper: safe logical -> FALSE
lf <- function(x) dplyr::coalesce(as.logical(x), FALSE)

# Define step flags explicitly (KEEP THESE STABLE ACROSS SITES)
consort_steps <- flags %>%
  mutate(
    step_01_candidates = TRUE,                                    # ICU candidates already in `flags`
    step_02_adult      = lf(adult),
    step_03_icu24h     = lf(icu_24h),
    step_04_has_county = !is.na(county_fips),                     # county_fips you created
    step_05_eligible   = lf(meets_data_rule),                     # ABG OR cont SpO2
    step_06_arf        = lf(arf_criterion_met)                    # phys ARF criteria
  ) %>%
  summarise(
    n_step_01_candidates = n(),
    n_step_02_adult      = sum(step_01_candidates & step_02_adult),
    n_step_03_icu24h     = sum(step_01_candidates & step_02_adult & step_03_icu24h),
    n_step_04_has_county = sum(step_01_candidates & step_02_adult & step_03_icu24h & step_04_has_county),
    n_step_05_eligible   = sum(step_01_candidates & step_02_adult & step_03_icu24h & step_04_has_county & step_05_eligible),
    n_step_06_arf        = sum(step_01_candidates & step_02_adult & step_03_icu24h & step_04_has_county & step_05_eligible & step_06_arf)
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "step", values_to = "remaining") %>%
  mutate(
    step_order = dplyr::case_when(
      step == "n_step_01_candidates" ~ 1L,
      step == "n_step_02_adult"      ~ 2L,
      step == "n_step_03_icu24h"     ~ 3L,
      step == "n_step_04_has_county" ~ 4L,
      step == "n_step_05_eligible"   ~ 5L,
      step == "n_step_06_arf"        ~ 6L,
      TRUE ~ NA_integer_
    ),
    step_label = dplyr::case_when(
      step == "n_step_01_candidates" ~ "ICU candidates (in years window)",
      step == "n_step_02_adult"      ~ "Age ≥ 18",
      step == "n_step_03_icu24h"     ~ "ICU LOS ≥ 24h",
      step == "n_step_04_has_county" ~ "County FIPS present",
      step == "n_step_05_eligible"   ~ "Eligible physiology (ABG or cont SpO2)",
      step == "n_step_06_arf"        ~ "Meets physiologic ARF criteria",
      TRUE ~ step
    )
  ) %>%
  arrange(step_order) %>%
  mutate(
    excluded_at_step = dplyr::lag(remaining, default = remaining[1]) - remaining,
    excluded_at_step = dplyr::if_else(step_order == 1L, NA_real_, as.numeric(excluded_at_step))
  ) %>%
  mutate(
    site_id = site_name,
    run_date = format(Sys.Date(), "%Y-%m-%d")
  ) %>%
  dplyr::select(site_id, run_date, step_order, step_label, remaining, excluded_at_step)

# Save (use your existing save_tbl() helper if present)
if (exists("save_tbl")) {
  save_tbl(consort_steps, "consort/consort_steps")
} else {
  out_dir2 <- file.path("output", "final", "consort")
  dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(consort_steps, file.path(out_dir2, paste0("consort_steps_", site_name, "_", format(Sys.Date(), "%Y%m%d"), ".csv")))
}

# =========================
# CONSORT exclusion reasons
# =========================

consort_reasons <- flags %>%
  mutate(
    fail_under18        = !lf(adult),
    fail_icu_lt24h       = lf(adult) & !lf(icu_24h),
    fail_missing_county  = lf(adult) & lf(icu_24h) & is.na(county_fips),
    fail_not_eligible    = lf(adult) & lf(icu_24h) & !is.na(county_fips) & !lf(meets_data_rule),
    fail_not_arf         = lf(adult) & lf(icu_24h) & !is.na(county_fips) & lf(meets_data_rule) & !lf(arf_criterion_met),
    
    # first-fail in step order (mutually exclusive)
    exclusion_reason = dplyr::case_when(
      fail_under18       ~ "Under 18",
      fail_icu_lt24h      ~ "ICU LOS < 24h",
      fail_missing_county ~ "Missing county FIPS",
      fail_not_eligible   ~ "No ABG or continuous SpO2 (±24h)",
      fail_not_arf        ~ "No physiologic ARF criteria (±24h)",
      TRUE                ~ "Included"
    )
  ) %>%
  count(exclusion_reason, name = "n") %>%
  mutate(
    site_id = site_name,
    run_date = format(Sys.Date(), "%Y-%m-%d")
  ) %>%
  arrange(desc(n))

if (exists("save_tbl")) {
  save_tbl(consort_reasons, "consort/consort_reasons")
} else {
  out_dir2 <- file.path("output", "final", "consort")
  dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(consort_reasons, file.path(out_dir2, paste0("consort_reasons_", site_name, "_", format(Sys.Date(), "%Y%m%d"), ".csv")))
}

# =========================
# TABLE 1: ARF cohort descriptive statistics
# =========================

out_dir_table1 <- file.path("output", "final", "table1")
dir.create(out_dir_table1, recursive = TRUE, showWarnings = FALSE)

# ARF analytic cohort (1 row per hospitalization)
table1_dat <- flags %>%
  filter(include_arf) %>%
  mutate(
    sex_category = coalesce(as.character(sex_category), "Missing"),
    race_category = coalesce(as.character(race_category), "Missing"),
    ethnicity_category = coalesce(as.character(ethnicity_category), "Missing"),
    age_band = coalesce(as.character(age_band), "Missing"),
    county_fips = coalesce(as.character(county_fips), "Missing")
  )

# -------------------------
# 1) Overall summary
# -------------------------
table1_overall <- tibble(
  site_name = site_name,
  n_arf = nrow(table1_dat),
  
  age_mean = mean(table1_dat$age_years, na.rm = TRUE),
  age_sd = sd(table1_dat$age_years, na.rm = TRUE),
  age_median = median(table1_dat$age_years, na.rm = TRUE),
  age_p25 = unname(quantile(table1_dat$age_years, 0.25, na.rm = TRUE)),
  age_p75 = unname(quantile(table1_dat$age_years, 0.75, na.rm = TRUE)),
  
  icu_los_hours_mean = mean(table1_dat$icu_los_hours, na.rm = TRUE),
  icu_los_hours_sd = sd(table1_dat$icu_los_hours, na.rm = TRUE),
  icu_los_hours_median = median(table1_dat$icu_los_hours, na.rm = TRUE),
  icu_los_hours_p25 = unname(quantile(table1_dat$icu_los_hours, 0.25, na.rm = TRUE)),
  icu_los_hours_p75 = unname(quantile(table1_dat$icu_los_hours, 0.75, na.rm = TRUE)),
  
  n_hypoxemic = sum(table1_dat$any_hypox, na.rm = TRUE),
  pct_hypoxemic = 100 * mean(table1_dat$any_hypox, na.rm = TRUE),
  
  n_hypercapnic = sum(table1_dat$any_hypercap, na.rm = TRUE),
  pct_hypercapnic = 100 * mean(table1_dat$any_hypercap, na.rm = TRUE),
  
  n_mixed = sum(table1_dat$mixed_arf, na.rm = TRUE),
  pct_mixed = 100 * mean(table1_dat$mixed_arf, na.rm = TRUE),
  
  n_in_hosp_death = sum(table1_dat$in_hosp_death, na.rm = TRUE),
  pct_in_hosp_death = 100 * mean(table1_dat$in_hosp_death, na.rm = TRUE),
  
  n_30d_death = sum(table1_dat$death_30d, na.rm = TRUE),
  pct_30d_death = 100 * mean(table1_dat$death_30d, na.rm = TRUE),
  
  n_has_abg = sum(table1_dat$has_abg, na.rm = TRUE),
  pct_has_abg = 100 * mean(table1_dat$has_abg, na.rm = TRUE),
  
  n_has_cont_spo2 = sum(table1_dat$has_cont_spo2, na.rm = TRUE),
  pct_has_cont_spo2 = 100 * mean(table1_dat$has_cont_spo2, na.rm = TRUE)
)

readr::write_csv(
  table1_overall,
  file.path(out_dir_table1, glue("table1_overall_{site_name}_{export_stamp}.csv")),
  na = ""
)

# -------------------------
# 2) Categorical distributions
# -------------------------
make_cat_tbl <- function(df, var, var_name) {
  df %>%
    count(value = .data[[var]], name = "n") %>%
    mutate(
      pct = 100 * n / sum(n),
      variable = var_name,
      site_name = site_name
    ) %>%
    select(site_name, variable, value, n, pct) %>%
    arrange(desc(n))
}

table1_sex <- make_cat_tbl(table1_dat, "sex_category", "sex_category")
table1_race <- make_cat_tbl(table1_dat, "race_category", "race_category")
table1_ethnicity <- make_cat_tbl(table1_dat, "ethnicity_category", "ethnicity_category")
table1_age_band <- make_cat_tbl(table1_dat, "age_band", "age_band")

readr::write_csv(
  table1_sex,
  file.path(out_dir_table1, glue("table1_sex_{site_name}_{export_stamp}.csv")),
  na = ""
)

readr::write_csv(
  table1_race,
  file.path(out_dir_table1, glue("table1_race_{site_name}_{export_stamp}.csv")),
  na = ""
)

readr::write_csv(
  table1_ethnicity,
  file.path(out_dir_table1, glue("table1_ethnicity_{site_name}_{export_stamp}.csv")),
  na = ""
)

readr::write_csv(
  table1_age_band,
  file.path(out_dir_table1, glue("table1_age_band_{site_name}_{export_stamp}.csv")),
  na = ""
)

# -------------------------
# 3) Formatted long Table 1
# -------------------------
fmt_num <- function(x, digits = 1) format(round(x, digits), nsmall = digits, trim = TRUE)

table1_long <- bind_rows(
  tibble(
    site_name = site_name,
    characteristic = "ARF hospitalizations, n",
    value = as.character(nrow(table1_dat))
  ),
  tibble(
    site_name = site_name,
    characteristic = "Age, mean (SD), y",
    value = paste0(fmt_num(mean(table1_dat$age_years, na.rm = TRUE)),
                   " (", fmt_num(sd(table1_dat$age_years, na.rm = TRUE)), ")")
  ),
  tibble(
    site_name = site_name,
    characteristic = "Age, median [IQR], y",
    value = paste0(
      fmt_num(median(table1_dat$age_years, na.rm = TRUE)),
      " [",
      fmt_num(quantile(table1_dat$age_years, 0.25, na.rm = TRUE)),
      ", ",
      fmt_num(quantile(table1_dat$age_years, 0.75, na.rm = TRUE)),
      "]"
    )
  ),
  tibble(
    site_name = site_name,
    characteristic = "ICU LOS, median [IQR], h",
    value = paste0(
      fmt_num(median(table1_dat$icu_los_hours, na.rm = TRUE)),
      " [",
      fmt_num(quantile(table1_dat$icu_los_hours, 0.25, na.rm = TRUE)),
      ", ",
      fmt_num(quantile(table1_dat$icu_los_hours, 0.75, na.rm = TRUE)),
      "]"
    )
  ),
  tibble(
    site_name = site_name,
    characteristic = "Hypoxemic ARF, n (%)",
    value = paste0(
      sum(table1_dat$any_hypox, na.rm = TRUE), " (",
      fmt_num(100 * mean(table1_dat$any_hypox, na.rm = TRUE)), "%)"
    )
  ),
  tibble(
    site_name = site_name,
    characteristic = "Hypercapnic ARF, n (%)",
    value = paste0(
      sum(table1_dat$any_hypercap, na.rm = TRUE), " (",
      fmt_num(100 * mean(table1_dat$any_hypercap, na.rm = TRUE)), "%)"
    )
  ),
  tibble(
    site_name = site_name,
    characteristic = "Mixed ARF, n (%)",
    value = paste0(
      sum(table1_dat$mixed_arf, na.rm = TRUE), " (",
      fmt_num(100 * mean(table1_dat$mixed_arf, na.rm = TRUE)), "%)"
    )
  ),
  tibble(
    site_name = site_name,
    characteristic = "In-hospital death, n (%)",
    value = paste0(
      sum(table1_dat$in_hosp_death, na.rm = TRUE), " (",
      fmt_num(100 * mean(table1_dat$in_hosp_death, na.rm = TRUE)), "%)"
    )
  ),
  tibble(
    site_name = site_name,
    characteristic = "30-day death, n (%)",
    value = paste0(
      sum(table1_dat$death_30d, na.rm = TRUE), " (",
      fmt_num(100 * mean(table1_dat$death_30d, na.rm = TRUE)), "%)"
    )
  )
)

readr::write_csv(
  table1_long,
  file.path(out_dir_table1, glue("table1_formatted_{site_name}_{export_stamp}.csv")),
  na = ""
)

message("✅ Wrote Table 1 outputs:")
message("  - ", file.path(out_dir_table1, glue("table1_overall_{site_name}_{export_stamp}.csv")))
message("  - ", file.path(out_dir_table1, glue("table1_sex_{site_name}_{export_stamp}.csv")))
message("  - ", file.path(out_dir_table1, glue("table1_race_{site_name}_{export_stamp}.csv")))
message("  - ", file.path(out_dir_table1, glue("table1_ethnicity_{site_name}_{export_stamp}.csv")))
message("  - ", file.path(out_dir_table1, glue("table1_age_band_{site_name}_{export_stamp}.csv")))
message("  - ", file.path(out_dir_table1, glue("table1_formatted_{site_name}_{export_stamp}.csv")))


message("\n",
        "------------------------------------------------------------------\n",
        "✔ CLIF site data extraction complete: ", site_name, "\n",
        "------------------------------------------------------------------\n",
        "\n",
        "All aggregated county-level outputs have been successfully generated.\n",
        "\n",
        "Action required:\n",
        "Upload the full 'output/final' directory to the shared Box folder.\n",
        "\n",
        "Only aggregated data are included. No patient-level information is present.\n",
        "\n",
        "If the script terminated without errors, your site’s contribution is complete.\n",
        "------------------------------------------------------------------\n"
)
