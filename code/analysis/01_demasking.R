# ================================================================================================
# COORDINATOR DEMASKING CHUNK
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

mask_version <- "mask_v20260320"
key_root <- file.path("offset_keys", mask_version, "master")
returned_dir <- file.path("returned_site_files")

master_key_county_year <- readr::read_csv(
  file.path(key_root, paste0("master_offset_site_county_year_", mask_version, ".csv")),
  show_col_types = FALSE
)

master_key_county_year_age_sex <- readr::read_csv(
  file.path(key_root, paste0("master_offset_site_county_year_age_sex_", mask_version, ".csv")),
  show_col_types = FALSE
)

# -------------------------------
# Pool and demask county-year
# -------------------------------

county_year_files <- list.files(
  returned_dir,
  pattern = "^site_county_year_.*\\.csv$",
  full.names = TRUE
)

pooled_county_year <- purrr::map_dfr(county_year_files, readr::read_csv, show_col_types = FALSE) %>%
  group_by(county_fips, year, mask_version) %>%
  summarise(
    A_all_ct  = sum(A_all_ct),
    A_elig_ct = sum(A_elig_ct),
    Y_ct      = sum(Y_ct),
    D_ct      = sum(D_ct),
    D30_ct    = sum(D30_ct),
    A_all_missing_county  = sum(A_all_missing_county),
    A_elig_missing_county = sum(A_elig_missing_county),
    Y_missing_county      = sum(Y_missing_county),
    D_missing_county      = sum(D_missing_county),
    .groups = "drop"
  ) %>%
  left_join(
    master_key_county_year %>%
      select(county_fips, year, mask_version, master_offset_n),
    by = c("county_fips", "year", "mask_version")
  ) %>%
  mutate(
    A_all_ct              = A_all_ct - master_offset_n,
    A_elig_ct             = A_elig_ct - master_offset_n,
    Y_ct                  = Y_ct - master_offset_n,
    D_ct                  = D_ct - master_offset_n,
    D30_ct                = D30_ct - master_offset_n,
    A_all_missing_county  = A_all_missing_county - master_offset_n,
    A_elig_missing_county = A_elig_missing_county - master_offset_n,
    Y_missing_county      = Y_missing_county - master_offset_n,
    D_missing_county      = D_missing_county - master_offset_n
  )

# -------------------------------
# Pool and demask county-year-age-sex
# -------------------------------

county_year_age_sex_files <- list.files(
  returned_dir,
  pattern = "^site_county_year_age_sex_.*\\.csv$",
  full.names = TRUE
)

pooled_county_year_age_sex <- purrr::map_dfr(county_year_age_sex_files, readr::read_csv, show_col_types = FALSE) %>%
  group_by(county_fips, year, age_band, sex, mask_version) %>%
  summarise(
    A_all_ctg  = sum(A_all_ctg),
    A_elig_ctg = sum(A_elig_ctg),
    Y_ctg      = sum(Y_ctg),
    D_ctg      = sum(D_ctg),
    .groups = "drop"
  ) %>%
  left_join(
    master_key_county_year_age_sex %>%
      select(county_fips, year, age_band, sex, mask_version, master_offset_n),
    by = c("county_fips", "year", "age_band", "sex", "mask_version")
  ) %>%
  mutate(
    A_all_ctg  = A_all_ctg - master_offset_n,
    A_elig_ctg = A_elig_ctg - master_offset_n,
    Y_ctg      = Y_ctg - master_offset_n,
    D_ctg      = D_ctg - master_offset_n
  )