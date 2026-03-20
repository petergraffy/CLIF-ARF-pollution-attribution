# ================================================================================================
# CLIF ARF National Estimation | Offset Key Generator (Coordinator Only)
# Creates:
#   1) site-specific keys for site_county_year
#   2) site-specific keys for site_county_year_age_sex
#   3) master keys (sum of all 10 site-specific offsets) for each table
#
# IMPORTANT:
# - With min_site_offset = 11 and n_sites = 10, master offsets are 110-150
# - These keys should NEVER be distributed except the site-specific file for each site
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# ------------------------------- Config ----------------------------------------------------------

mask_version <- "mask_v20260320"
site_names <- c(
  "site01", "site02", "site03", "site04", "site05",
  "site06", "site07", "site08", "site09", "site10"
)

years <- 2018:2024
age_bands <- c("18-39", "40-64", "65-74", "75+")
sexes <- c("Female", "Male")

n_sites <- length(site_names)
min_site_offset <- 11L   # strictly >10
extra_total_max <- 40L   # additional total per cell distributed across sites

# Directories
key_root <- file.path("offset_keys", mask_version)
dir.create(key_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(key_root, "site_county_year"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(key_root, "site_county_year_age_sex"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(key_root, "master"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------- Helpers ---------------------------------------------------------

# CONUS = drop AK, HI, and territories
# This uses the built-in 'maps' package county crosswalk, which is lightweight and usually available.
get_conus_counties <- function() {
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop(
      "Package 'maps' is required for county FIPS generation.\n",
      "Install once with: install.packages('maps')"
    )
  }
  
  county_fips <- maps::county.fips %>%
    as_tibble() %>%
    transmute(
      county_fips = sprintf("%05d", fips)
    ) %>%
    distinct() %>%
    filter(
      !substr(county_fips, 1, 2) %in% c(
        "02", # Alaska
        "15", # Hawaii
        "60", # American Samoa
        "66", # Guam
        "69", # Northern Mariana Islands
        "72", # Puerto Rico
        "78"  # U.S. Virgin Islands
      )
    ) %>%
    arrange(county_fips)
  
  county_fips
}

# Integer partition of 'extra_total' across n_sites
# Returns integer vector length n_sites summing exactly to extra_total
partition_extra <- function(extra_total, n_sites) {
  if (extra_total == 0L) return(rep.int(0L, n_sites))
  as.integer(rmultinom(n = 1, size = extra_total, prob = rep(1, n_sites)))
}

# Generate site-specific offsets for one scaffold
# One offset per cell per site; same offset will be added to all masked count columns for that table
generate_site_offsets <- function(scaffold_df, site_names, min_site_offset = 11L,
                                  extra_total_max = 40L, seed = 1L) {
  set.seed(seed)
  
  n_cells <- nrow(scaffold_df)
  n_sites <- length(site_names)
  
  extra_totals <- sample.int(extra_total_max + 1L, size = n_cells, replace = TRUE) - 1L
  
  extras_mat <- matrix(0L, nrow = n_cells, ncol = n_sites)
  for (i in seq_len(n_cells)) {
    extras_mat[i, ] <- partition_extra(extra_totals[i], n_sites)
  }
  
  offsets_mat <- extras_mat + min_site_offset
  
  colnames(offsets_mat) <- site_names
  
  offsets_long <- as_tibble(offsets_mat) %>%
    mutate(.row_id = seq_len(n_cells)) %>%
    pivot_longer(
      cols = all_of(site_names),
      names_to = "site_name",
      values_to = "offset_n"
    )
  
  scaffold_df %>%
    mutate(.row_id = seq_len(n())) %>%
    left_join(offsets_long, by = ".row_id") %>%
    mutate(
      offset_n = as.integer(offset_n),
      mask_version = mask_version
    ) %>%
    select(-.row_id)
}

# ------------------------------- Build scaffolds -------------------------------------------------

conus_counties <- get_conus_counties()

# Table 1 scaffold: county_fips x year
scaffold_county_year <- tidyr::expand_grid(
  county_fips = conus_counties$county_fips,
  year = years
) %>%
  mutate(
    cell_id = paste(county_fips, year, sep = "__")
  ) %>%
  select(cell_id, county_fips, year)

# Table 2 scaffold: county_fips x year x age_band x sex
scaffold_county_year_age_sex <- tidyr::expand_grid(
  county_fips = conus_counties$county_fips,
  year = years,
  age_band = age_bands,
  sex = sexes
) %>%
  mutate(
    cell_id = paste(county_fips, year, age_band, sex, sep = "__")
  ) %>%
  select(cell_id, county_fips, year, age_band, sex)

message("Scaffold sizes:")
message("  county_year: ", nrow(scaffold_county_year))
message("  county_year_age_sex: ", nrow(scaffold_county_year_age_sex))

# ------------------------------- Generate site-specific keys -------------------------------------

site_keys_county_year <- generate_site_offsets(
  scaffold_df = scaffold_county_year,
  site_names = site_names,
  min_site_offset = min_site_offset,
  extra_total_max = extra_total_max,
  seed = 1001L
)

site_keys_county_year_age_sex <- generate_site_offsets(
  scaffold_df = scaffold_county_year_age_sex,
  site_names = site_names,
  min_site_offset = min_site_offset,
  extra_total_max = extra_total_max,
  seed = 2001L
)

# ------------------------------- Master keys -----------------------------------------------------

master_key_county_year <- site_keys_county_year %>%
  group_by(cell_id, county_fips, year, mask_version) %>%
  summarise(
    master_offset_n = sum(offset_n),
    .groups = "drop"
  ) %>%
  arrange(year, county_fips)

master_key_county_year_age_sex <- site_keys_county_year_age_sex %>%
  group_by(cell_id, county_fips, year, age_band, sex, mask_version) %>%
  summarise(
    master_offset_n = sum(offset_n),
    .groups = "drop"
  ) %>%
  arrange(year, county_fips, age_band, sex)

# ------------------------------- Write site-specific keys ----------------------------------------

walk(site_names, function(s) {
  key1 <- site_keys_county_year %>%
    filter(site_name == s) %>%
    arrange(year, county_fips)
  
  key2 <- site_keys_county_year_age_sex %>%
    filter(site_name == s) %>%
    arrange(year, county_fips, age_band, sex)
  
  write_csv(
    key1,
    file.path(key_root, "site_county_year", paste0("offset_key_site_county_year_", s, "_", mask_version, ".csv")),
    na = ""
  )
  
  write_csv(
    key2,
    file.path(key_root, "site_county_year_age_sex", paste0("offset_key_site_county_year_age_sex_", s, "_", mask_version, ".csv")),
    na = ""
  )
})

# ------------------------------- Write master keys ------------------------------------------------

write_csv(
  master_key_county_year,
  file.path(key_root, "master", paste0("master_offset_site_county_year_", mask_version, ".csv")),
  na = ""
)

write_csv(
  master_key_county_year_age_sex,
  file.path(key_root, "master", paste0("master_offset_site_county_year_age_sex_", mask_version, ".csv")),
  na = ""
)

# ------------------------------- Write audit/QC summary ------------------------------------------

audit_summary <- bind_rows(
  site_keys_county_year %>%
    group_by(site_name, mask_version) %>%
    summarise(
      table_name = "site_county_year",
      n_rows = n(),
      min_offset = min(offset_n),
      max_offset = max(offset_n),
      sum_offset = sum(offset_n),
      .groups = "drop"
    ),
  site_keys_county_year_age_sex %>%
    group_by(site_name, mask_version) %>%
    summarise(
      table_name = "site_county_year_age_sex",
      n_rows = n(),
      min_offset = min(offset_n),
      max_offset = max(offset_n),
      sum_offset = sum(offset_n),
      .groups = "drop"
    )
)

master_summary <- bind_rows(
  master_key_county_year %>%
    summarise(
      table_name = "site_county_year",
      n_rows = n(),
      min_master_offset = min(master_offset_n),
      max_master_offset = max(master_offset_n),
      sum_master_offset = sum(master_offset_n)
    ),
  master_key_county_year_age_sex %>%
    summarise(
      table_name = "site_county_year_age_sex",
      n_rows = n(),
      min_master_offset = min(master_offset_n),
      max_master_offset = max(master_offset_n),
      sum_master_offset = sum(master_offset_n)
    )
)

write_csv(
  audit_summary,
  file.path(key_root, "master", paste0("site_key_audit_", mask_version, ".csv")),
  na = ""
)

write_csv(
  master_summary,
  file.path(key_root, "master", paste0("master_key_summary_", mask_version, ".csv")),
  na = ""
)

message("✅ Offset keys created successfully")
message("Key root: ", normalizePath(key_root))
message("Expected master offset range with min_site_offset=11 and extra_total_max=40: 110 to 150")




