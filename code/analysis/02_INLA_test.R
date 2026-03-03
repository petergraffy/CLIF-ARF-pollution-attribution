# =====================================================================================
# Illinois dry run: UChicago export -> INLA-SPDE ARF per-capita model + g-computation
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(tidyr)
  library(sf)
  library(tigris)
  library(tidycensus)
  library(INLA)       # make sure INLA is installed
  library(ggplot2)
})

options(tigris_use_cache = TRUE)
sf_use_s2(FALSE)

# ---------------------------
# 0) USER INPUTS
# ---------------------------

# Path to your site exports (UChicago)
# (use the newest timestamped files you just generated)
path_cty_year <- "output/site_county_year_UCMC_20260303_092354.csv"
path_cty_year_as <- "output/site_county_year_age_sex_UCMC_20260303_092354.csv" # optional

# Your county-year pollution files (you already have these in your exposome folder)
# Expected columns: county_fips (5-char), year (int), pm25, no2 (or similar)
path_pm25 <- "exposome/pm25_county_year.csv"
path_no2  <- "exposome/no2_county_year.csv"

# Census API key (set once per machine). If already set, you can skip.
# census_api_key("YOUR_KEY_HERE", install = TRUE, overwrite = FALSE)

years_use <- 2018:2024   # set to 2018:2025 once your export includes 2025

# ---------------------------
# 1) READ EXPORTS
# ---------------------------

dat <- read_csv(path_cty_year, show_col_types = FALSE) %>%
  mutate(
    county_fips = str_pad(as.character(county_fips), 5, pad="0"),
    year = as.integer(year)
  )

# Sanity: required columns (after your “both denominators” update)
stopifnot(all(c("A_all_ct","A_elig_ct","Y_ct","D_ct","year","county_fips") %in% names(dat)))

# IL counties (FIPS state=17)
dat_il <- dat %>%
  filter(str_sub(county_fips, 1, 2) == "17") %>%
  filter(year %in% years_use)

stopifnot(nrow(dat_il) > 0)

# Phenotype observability omega (avoid division by 0)
dat_il <- dat_il %>%
  mutate(
    omega = ifelse(A_all_ct > 0, A_elig_ct / A_all_ct, NA_real_),
    omega = pmin(pmax(omega, 1e-6), 1)  # clamp for log-offset stability
  )

# ---------------------------
# 2) ADD ACS POPULATION (N_ct)
# ---------------------------

# Total population per county-year (ACS 5-year is fine for test)
# Note: ACS "year" corresponds to release year. We'll request each year.
get_acs_pop <- function(y) {
  tidycensus::get_acs(
    geography = "county",
    variables = c(pop = "B01003_001"),
    year = y,
    survey = "acs5",
    cache_table = TRUE
  ) %>%
    transmute(
      county_fips = GEOID,
      year = y,
      N = estimate
    )
}

acs_pop <- bind_rows(lapply(years_use, get_acs_pop)) %>%
  filter(str_sub(county_fips, 1, 2) == "17")

dat_il <- dat_il %>%
  left_join(acs_pop, by = c("county_fips","year"))

stopifnot(sum(is.na(dat_il$N)) == 0)

# ---------------------------
# 3) ADD POLLUTION (PM2.5, NO2)
# ---------------------------

pm25 <- read_csv(path_pm25, show_col_types = FALSE) %>%
  transmute(county_fips = str_pad(as.character(GEOID), 5, pad="0"),
            year = as.integer(year),
            pm25 = as.numeric(pm25_mean))

no2 <- read_csv(path_no2, show_col_types = FALSE) %>%
  transmute(county_fips = str_pad(as.character(GEOID), 5, pad="0"),
            year = as.integer(year),
            no2 = as.numeric(no2_mean))

dat_il <- dat_il %>%
  left_join(pm25, by = c("county_fips","year")) %>%
  left_join(no2,  by = c("county_fips","year"))

stopifnot(sum(is.na(dat_il$pm25)) == 0, sum(is.na(dat_il$no2)) == 0)

# Counterfactual pollution levels (10th percentile within IL-years used; you can switch to national later)
pm25_cf <- unname(quantile(dat_il$pm25, probs = 0.10, na.rm = TRUE))
no2_cf  <- unname(quantile(dat_il$no2,  probs = 0.10, na.rm = TRUE))

# ---------------------------
# 4) GEOMETRY + IL-ONLY MESH
# ---------------------------

il_counties <- tigris::counties(state = "IL", cb = TRUE, year = 2022) %>%
  st_as_sf() %>%
  transmute(county_fips = GEOID, geometry)

# Join geometry to data (for centroids + mapping)
dat_il_sf <- il_counties %>%
  left_join(dat_il, by = "county_fips") %>%
  filter(!is.na(year))  # keep only county-years present

# Use a projected CRS for mesh
dat_il_sf <- st_transform(dat_il_sf, 5070)  # NAD83 / Conus Albers

cent <- st_centroid(dat_il_sf) %>%
  st_coordinates()

# One centroid per row (county-year rows repeat: OK for A-matrix)
locs <- cent

# Build mesh over IL
mesh <- inla.mesh.2d(
  loc = locs,
  max.edge = c(50e3, 150e3),   # tune if needed
  cutoff   = 10e3
)

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(100e3, 0.5), # P(range < 100km) = 0.5
  prior.sigma = c(1, 0.01)     # P(sigma > 1) = 0.01
)

A_mat <- inla.spde.make.A(mesh, loc = locs)
s_index <- inla.spde.make.index("s", spde$n.spde)

# ---------------------------
# 5) MODEL MATRIX + OFFSETS
# ---------------------------

# Outcome: Y_ct (ARF cases)
# Offset: log(N_ct) + log(omega_ct)
dat_model <- dat_il_sf %>%
  st_drop_geometry() %>%
  mutate(
    year_f = as.integer(year),
    offset_log = log(N) + log(omega)
  )

# INLA stack
stk <- inla.stack(
  data = list(y = dat_model$Y_ct),
  A = list(A_mat, 1),
  effects = list(
    s = s_index,
    data.frame(
      intercept = 1,
      pm25 = scale(dat_model$pm25)[,1],
      no2  = scale(dat_model$no2)[,1],
      year_id = as.integer(factor(dat_model$year_f))
    )
  ),
  tag = "est"
)

# Formula: NegBin (overdispersion), SPDE spatial + RW1 year
formula <- y ~ 0 + intercept + pm25 + no2 +
  f(s, model = spde) +
  f(year_id, model = "rw1")

fit <- inla(
  formula,
  family = "nbinomial",
  data = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk), compute = TRUE, link = 1),
  control.fixed = list(expand.factor.strategy = "inla"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  offset = dat_model$offset_log
)

print(fit$summary.fixed)
cat("\nDIC:", fit$dic$dic, " WAIC:", fit$waic$waic, "\n")

# ---------------------------
# 6) COUNTERFACTUAL PREDICTIONS (g-computation)
# ---------------------------

# Build two new datasets: observed vs counterfactual pollution
dat_obs <- dat_model %>%
  mutate(pm25 = pm25, no2 = no2)

dat_cf <- dat_model %>%
  mutate(pm25 = pm25_cf, no2 = no2_cf)

# NOTE: we scaled pm25/no2 in the stack. Recompute scaling using original model scaling.
pm25_center <- mean(dat_model$pm25, na.rm = TRUE)
pm25_scale  <- sd(dat_model$pm25, na.rm = TRUE)
no2_center  <- mean(dat_model$no2,  na.rm = TRUE)
no2_scale   <- sd(dat_model$no2,  na.rm = TRUE)

make_lp <- function(pm25_val, no2_val) {
  (pm25_val - pm25_center)/pm25_scale
}

pm25_obs_sc <- (dat_obs$pm25 - pm25_center)/pm25_scale
no2_obs_sc  <- (dat_obs$no2  - no2_center)/no2_scale
pm25_cf_sc  <- (dat_cf$pm25  - pm25_center)/pm25_scale
no2_cf_sc   <- (dat_cf$no2   - no2_center)/no2_scale

# Posterior linear predictor for observed and CF using fixed effects only is insufficient because spatial+year matter.
# Instead: use fitted predictor means and then adjust by delta from fixed effects (approx).
# For a dry run, do a simple approach: predict expected counts using fitted marginals at observed,
# then shift by beta*(x_cf - x_obs) on the log scale.

beta_pm25 <- fit$summary.fixed["pm25","mean"]
beta_no2  <- fit$summary.fixed["no2","mean"]

eta_shift <- beta_pm25*(pm25_cf_sc - pm25_obs_sc) + beta_no2*(no2_cf_sc - no2_obs_sc)

# indices for the estimation stack rows
idx_est <- inla.stack.index(stk, tag = "est")$data

mu_obs <- fit$summary.fitted.values$mean[idx_est]

stopifnot(length(mu_obs) == nrow(dat_model))  # should now be TRUE

mu_cf <- mu_obs * exp(eta_shift)

dat_out <- dat_model %>%
  mutate(
    mu_obs = mu_obs,
    mu_cf  = mu_cf,
    excess = pmax(mu_obs - mu_cf, 0)
  )

# Aggregate to IL totals by year
il_year <- dat_out %>%
  group_by(year) %>%
  summarise(
    arf_obs = sum(mu_obs, na.rm = TRUE),
    arf_cf  = sum(mu_cf,  na.rm = TRUE),
    excess  = sum(excess, na.rm = TRUE),
    .groups = "drop"
  )

print(il_year)

# ---------------------------
# 7) MAPS (quick)
# ---------------------------

map_df <- il_counties %>%
  left_join(
    dat_out %>% group_by(county_fips) %>%
      summarise(excess = sum(excess, na.rm = TRUE),
                arf_obs = sum(mu_obs, na.rm = TRUE),
                .groups="drop"),
    by = "county_fips"
  ) %>%
  st_as_sf()

map_df <- il_counties %>%
  left_join(
    dat_out %>%
      group_by(county_fips) %>%
      summarise(
        excess  = sum(excess, na.rm = TRUE),
        arf_obs = sum(mu_obs, na.rm = TRUE),
        pop     = sum(N, na.rm = TRUE),   # <-- add this
        .groups = "drop"
      ) %>%
      mutate(
        excess_per_100k = 1e5 * excess / pop,
        arf_rate_per_100k = 1e5 * arf_obs / pop
      ),
    by = "county_fips"
  ) %>%
  st_as_sf()

map_df <- il_counties %>%
  left_join(
    dat_out %>%
      group_by(county_fips) %>%
      summarise(
        excess  = sum(excess, na.rm = TRUE),
        arf_obs = sum(mu_obs, na.rm = TRUE),
        pop_mean = mean(N, na.rm = TRUE),
        n_years  = n_distinct(year),
        .groups = "drop"
      ) %>%
      mutate(
        excess_per_100k_annual = 1e5 * (excess / n_years) / pop_mean
      ),
    by = "county_fips"
  ) %>%
  st_as_sf()

p1 <- ggplot(map_df) +
  geom_sf(aes(fill = excess_per_100k_annual), color = "white", linewidth = 0.1) +
  labs(title = "Illinois: Estimated excess ARF cases attributable to pollution",
       subtitle = "UChicago-only dry run (pipeline validation)",
       fill = "Excess ARF per 100k") +
  theme_minimal(base_size = 12)

print(p1)

out_dir <- "output/il_test_run"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = file.path(out_dir, "il_excess_arf_map.png"),
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)



# Illinois boundary (projected)
il_boundary <- st_union(st_transform(il_counties, 5070))

# Regular grid (tune cellsize: 10km is a nice start)
grid <- st_make_grid(il_boundary, cellsize = 10e3, what = "centers") %>%
  st_as_sf() %>%
  st_intersection(il_boundary)

grid_xy <- st_coordinates(grid)

# Projector matrix from mesh nodes to grid points
A_grid <- inla.spde.make.A(mesh, loc = grid_xy)

# Extract posterior mean of the spatial random field at mesh nodes
s_mean <- fit$summary.random$s$mean  # length = spde$n.spde

# Interpolate to grid
grid$spatial_mean <- as.vector(A_grid %*% s_mean)

grid_df <- cbind(
  st_drop_geometry(grid),
  as.data.frame(grid_xy)
)
names(grid_df)[(ncol(grid_df)-1):ncol(grid_df)] <- c("x","y")

ggplot(grid_df, aes(x = x, y = y, fill = spatial_mean)) +
  geom_raster() +
  coord_equal() +
  labs(
    title = "Illinois: smooth spatial effect (INLA-SPDE)",
    subtitle = "Posterior mean of spatial random field",
    fill = "Spatial effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    panel.grid = element_blank()
  )

ggplot() +
  geom_raster(data = grid_df, aes(x = x, y = y, fill = spatial_mean)) +
  geom_sf(data = st_cast(il_boundary, "MULTIPOLYGON"), fill = NA, linewidth = 0.4) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    panel.grid = element_blank()
  )







