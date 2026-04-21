# =====================================================================================
# Illinois dry run (multi-site): pooled demasked exports -> INLA-SPDE ARF per-capita + g-comp
# Point this script to the pooled demasked county-year file created after external unmasking.
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
  library(INLA)
  library(ggplot2)
})

options(tigris_use_cache = TRUE)
sf_use_s2(FALSE)

# ---------------------------
# 0) USER INPUTS
# ---------------------------

years_use <- 2018:2024

# Pollution files (your formats)
path_pm25 <- "exposome/pm25_county_year.csv"  # columns: GEOID, year, pm25_mean
path_no2  <- "exposome/no2_county_year.csv"   # columns: GEOID, year, no2_mean

# Pooled demasked county-year file
path_pooled_county_year <- file.path(
  "output",
  "pooled_demasked",
  "pooled_demasked_site_county_year_mask_v20260320.csv"
)

# Output folder for this dry run
out_dir <- file.path("output", "il_test_run_2sites")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------
# 1) READ DEMASKED POOLED EXPORT
# ---------------------------

if (!file.exists(path_pooled_county_year)) {
  stop("Missing pooled demasked file: ", path_pooled_county_year,
       ". Create or place the pooled demasked file first.")
}

dat_il <- read_csv(path_pooled_county_year, show_col_types = FALSE) %>%
  mutate(
    county_fips = str_pad(as.character(county_fips), 5, pad = "0"),
    year = as.integer(year)
  ) %>%
  filter(str_sub(county_fips, 1, 2) == "17") %>%
  filter(year %in% years_use) %>%
  mutate(
    omega = ifelse(A_all_ct > 0, A_elig_ct / A_all_ct, NA_real_),
    omega = pmin(pmax(omega, 1e-6), 1)
  )

req_cols <- c("A_all_ct","A_elig_ct","Y_ct","D_ct","year","county_fips")
stopifnot(all(req_cols %in% names(dat_il)))
stopifnot(nrow(dat_il) > 0)

# ---------------------------
# 2) ACS POPULATION (N_ct)
# ---------------------------

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
# 3) POLLUTION (PM2.5, NO2)
# ---------------------------

pm25 <- read_csv(path_pm25, show_col_types = FALSE) %>%
  transmute(
    county_fips = str_pad(as.character(GEOID), 5, pad="0"),
    year = as.integer(year),
    pm25 = as.numeric(pm25_mean)
  )

no2 <- read_csv(path_no2, show_col_types = FALSE) %>%
  transmute(
    county_fips = str_pad(as.character(GEOID), 5, pad="0"),
    year = as.integer(year),
    no2 = as.numeric(no2_mean)
  )

dat_il <- dat_il %>%
  left_join(pm25, by = c("county_fips","year")) %>%
  left_join(no2,  by = c("county_fips","year"))

stopifnot(sum(is.na(dat_il$pm25)) == 0, sum(is.na(dat_il$no2)) == 0)

pm25_cf <- unname(quantile(dat_il$pm25, probs = 0.10, na.rm = TRUE))
no2_cf  <- unname(quantile(dat_il$no2,  probs = 0.10, na.rm = TRUE))

# ---------------------------
# 4) GEOMETRY + IL-ONLY MESH
# ---------------------------

il_counties <- tigris::counties(state = "IL", cb = TRUE, year = 2022) %>%
  st_as_sf() %>%
  transmute(county_fips = GEOID, geometry)

dat_il_sf <- il_counties %>%
  left_join(dat_il, by = "county_fips") %>%
  filter(!is.na(year))

dat_il_sf <- st_transform(dat_il_sf, 5070)

locs <- st_coordinates(st_centroid(dat_il_sf))

mesh <- fmesher::fm_mesh_2d_inla(
  loc = locs,
  max.edge = c(50e3, 150e3),
  cutoff = 10e3
)

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(100e3, 0.5),
  prior.sigma = c(1, 0.01)
)

A_mat   <- inla.spde.make.A(mesh, loc = locs)
s_index <- inla.spde.make.index("s", spde$n.spde)

# ---------------------------
# 5) MODEL MATRIX + OFFSETS
# ---------------------------

dat_model <- dat_il_sf %>%
  st_drop_geometry() %>%
  mutate(
    year_f = as.integer(year),
    offset_log = log(N) + log(omega)
  )

# Scale exposures once (use the same scaling for g-comp)
pm25_center <- mean(dat_model$pm25, na.rm = TRUE)
pm25_scale  <- sd(dat_model$pm25,  na.rm = TRUE)
no2_center  <- mean(dat_model$no2,  na.rm = TRUE)
no2_scale   <- sd(dat_model$no2,   na.rm = TRUE)

pm25_sc <- (dat_model$pm25 - pm25_center) / pm25_scale
no2_sc  <- (dat_model$no2  - no2_center)  / no2_scale

stk <- inla.stack(
  data = list(y = dat_model$Y_ct),
  A = list(A_mat, 1),
  effects = list(
    s = s_index,
    data.frame(
      intercept = 1,
      pm25 = pm25_sc,
      no2  = no2_sc,
      year_id = as.integer(factor(dat_model$year_f))
    )
  ),
  tag = "est"
)

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
# 6) COUNTERFACTUAL (approx g-comp)
# ---------------------------

beta_pm25 <- fit$summary.fixed["pm25","mean"]
beta_no2  <- fit$summary.fixed["no2","mean"]

pm25_cf_sc <- (pm25_cf - pm25_center) / pm25_scale
no2_cf_sc  <- (no2_cf  - no2_center)  / no2_scale

eta_shift <- beta_pm25*(pm25_cf_sc - pm25_sc) + beta_no2*(no2_cf_sc - no2_sc)

idx_est <- inla.stack.index(stk, tag = "est")$data
mu_obs <- fit$summary.fitted.values$mean[idx_est]
stopifnot(length(mu_obs) == nrow(dat_model))

mu_cf <- mu_obs * exp(eta_shift)

dat_out <- dat_model %>%
  mutate(
    mu_obs = mu_obs,
    mu_cf  = mu_cf,
    excess = pmax(mu_obs - mu_cf, 0)
  )

il_year <- dat_out %>%
  group_by(year) %>%
  summarise(
    arf_obs = sum(mu_obs, na.rm = TRUE),
    arf_cf  = sum(mu_cf,  na.rm = TRUE),
    excess  = sum(excess, na.rm = TRUE),
    .groups = "drop"
  )

print(il_year)
write_csv(il_year, file.path(out_dir, "il_totals_by_year_2sites.csv"))

# ---------------------------
# 7) MAPS (annual excess per 100k)
# ---------------------------

map_df <- il_counties %>%
  left_join(
    dat_out %>%
      group_by(county_fips) %>%
      summarise(
        excess  = sum(excess, na.rm = TRUE),
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
  labs(
    title = "Illinois: Estimated excess ARF attributable to pollution",
    subtitle = "Pooled demasked dry run",
    fill = "Excess ARF per 100k (annual)"
  ) +
  theme_minimal(base_size = 12)

print(p1)

ggsave(
  filename = file.path(out_dir, "il_excess_arf_map_2sites.png"),
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)

# ---------------------------
# 8) SMOOTH SPATIAL FIELD PLOT (posterior mean of spatial random effect)
# ---------------------------

il_boundary <- st_union(st_transform(il_counties, 5070))

grid <- st_make_grid(il_boundary, cellsize = 10e3, what = "centers") %>%
  st_as_sf() %>%
  st_intersection(il_boundary)

grid_xy <- st_coordinates(grid)

A_grid <- inla.spde.make.A(mesh, loc = grid_xy)
s_mean <- fit$summary.random$s$mean
grid$spatial_mean <- as.vector(A_grid %*% s_mean)

grid_df <- cbind(st_drop_geometry(grid), as.data.frame(grid_xy))
names(grid_df)[(ncol(grid_df)-1):ncol(grid_df)] <- c("x","y")

il_outline <- st_cast(il_boundary, "MULTILINESTRING") %>% st_as_sf()
il_xy <- st_coordinates(il_outline)
il_xy <- as.data.frame(il_xy)

p_spatial <- ggplot() +
  geom_raster(data = grid_df, aes(x = x, y = y, fill = spatial_mean)) +
  geom_path(data = il_xy, aes(x = X, y = Y, group = L1), linewidth = 0.4) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())

print(p_spatial)

ggsave(
  filename = file.path(out_dir, "il_spatial_field_2sites.png"),
  plot = p_spatial,
  width = 10,
  height = 8,
  dpi = 300
)

message("\n✔ Two-site Illinois INLA dry run complete. Outputs saved to: ", out_dir, "\n")





# ---- pull scaling from your model object ----
pm25_center <- mean(dat_model$pm25, na.rm = TRUE)
pm25_scale  <- sd(dat_model$pm25,  na.rm = TRUE)

# ---- choose exposure grid and reference ----
x_grid <- seq(min(dat_model$pm25, na.rm = TRUE),
              max(dat_model$pm25, na.rm = TRUE),
              length.out = 200)

x_ref  <- unname(quantile(dat_model$pm25, 0.10, na.rm = TRUE))  # reference = 10th pct

# difference in *scaled* exposure from reference
dx_sc <- (x_grid - x_ref) / pm25_scale

# ---- posterior of beta_pm25 ----
m_beta <- fit$marginals.fixed$pm25

# For each x, compute posterior of log(RR) = beta * dx_sc,
# then summarize RR = exp(logRR)
rr_summ <- lapply(dx_sc, function(d) {
  m_logrr <- INLA::inla.tmarginal(function(b) b * d, m_beta)
  m_rr    <- INLA::inla.tmarginal(exp, m_logrr)
  qs      <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), m_rr)
  c(lo = qs[1], mid = qs[2], hi = qs[3])
})

rr_df_pm25 <- data.frame(
  pm25 = x_grid,
  do.call(rbind, rr_summ)
) %>%
  mutate(ref = x_ref)

p_pm25 <- ggplot(rr_df_pm25, aes(x = pm25, y = mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(
    title = "Posterior exposure–response: PM2.5",
    subtitle = paste0("Rate ratio vs reference (PM2.5 = ", round(x_ref, 2), ")"),
    x = "PM2.5 (annual mean)",
    y = "Rate ratio (ARF incidence)"
  ) +
  theme_minimal(base_size = 12)

print(p_pm25)


ggsave(file.path(out_dir, "exposure_response_pm25.png"), p_pm25, width = 8, height = 5, dpi = 300)
ggsave(file.path(out_dir, "exposure_response_pm25.pdf"), p_pm25, width = 8, height = 5)


no2_center <- mean(dat_model$no2, na.rm = TRUE)
no2_scale  <- sd(dat_model$no2,  na.rm = TRUE)

x_grid <- seq(min(dat_model$no2, na.rm = TRUE),
              max(dat_model$no2, na.rm = TRUE),
              length.out = 200)

x_ref <- unname(quantile(dat_model$no2, 0.10, na.rm = TRUE))
dx_sc <- (x_grid - x_ref) / no2_scale

m_beta <- fit$marginals.fixed$no2

rr_summ <- lapply(dx_sc, function(d) {
  m_logrr <- INLA::inla.tmarginal(function(b) b * d, m_beta)
  m_rr    <- INLA::inla.tmarginal(exp, m_logrr)
  qs      <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), m_rr)
  c(lo = qs[1], mid = qs[2], hi = qs[3])
})

rr_df_no2 <- data.frame(
  no2 = x_grid,
  do.call(rbind, rr_summ)
)

p_no2 <- ggplot(rr_df_no2, aes(x = no2, y = mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(
    title = "Posterior exposure–response: NO₂",
    subtitle = paste0("Rate ratio vs reference (NO₂ = ", round(x_ref, 2), ")"),
    x = "NO₂ (annual mean)",
    y = "Rate ratio (ARF incidence)"
  ) +
  theme_minimal(base_size = 12)

print(p_no2)

ggsave(file.path(out_dir, "exposure_response_pm25.png"), p_no2, width = 8, height = 5, dpi = 300)
ggsave(file.path(out_dir, "exposure_response_pm25.pdf"), p_no2, width = 8, height = 5)

rr_both <- bind_rows(
  rr_df_pm25 %>% transmute(exposure = pm25, lo, mid, hi, pollutant = "PM2.5"),
  rr_df_no2  %>% transmute(exposure = no2,  lo, mid, hi, pollutant = "NO₂")
)

p_both <- ggplot(rr_both, aes(x = exposure, y = mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ pollutant, scales = "free_x") +
  labs(
    title = "Posterior exposure–response curves",
    x = "Annual mean concentration",
    y = "Rate ratio (vs 10th percentile reference)"
  ) +
  theme_minimal(base_size = 12)

print(p_both)
ggsave(file.path(out_dir, "exposure_response_curves.png"), p_both, width = 10, height = 4.5, dpi = 300)
