# ============================================================================
# INDIA'S SOLAR TRANSITION AND COAL UTILIZATION
# Reproducible R coding sample | Data: Ember, Jan 2019–Nov 2025
#
# Research question:
# Are increases in solar and wind generation shares associated with lower
# utilization of installed coal capacity across Indian states?
#
# Observational analysis: associations are NOT causal displacement estimates.
#
# RUN FROM THE PROJECT ROOT:
#   source("india_solar_coal_coding_sample.R")
#
# REQUIRED INPUT:
#   data/raw/india_monthly_full_release_long_format.csv
#
# OUTPUTS (created automatically):
#   data/processed/panel.rds, data/processed/regression_models.rds
#   figures/01_*.png through figures/06_*.png
#   tables/national_generation_mix.csv
#   tables/main_regression_results.csv
#   tables/leave_one_state_out.csv
#
# Dependencies: tidyverse, janitor, lubridate, fixest, fwildclusterboot, dqrng
# Install ordinary dependencies with install.packages(...).
# If fwildclusterboot is not available from CRAN, see the package's r-universe.
# ============================================================================

# 1. SETUP --------------------------------------------------------------------

required_packages <- c(
  "tidyverse", "janitor", "lubridate", "fixest",
  "fwildclusterboot", "dqrng"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the following R packages before running: ",
       paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

input_path <- file.path(
  "data", "raw", "india_monthly_full_release_long_format.csv"
)
if (!file.exists(input_path)) {
  stop("Input not found: ", input_path,
       "\nSet your working directory to the project root and place the Ember CSV in data/raw/.")
}

for (folder in c("data/processed", "figures", "tables")) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

save_figure <- function(plot, filename, width = 10, height = 7) {
  ggplot2::ggsave(
    filename = file.path("figures", filename),
    plot = plot, width = width, height = height, dpi = 300
  )
}

# 2. CLEAN AND HARMONIZE THE STATE-MONTH DATA -------------------------------

raw <- readr::read_csv(input_path, show_col_types = FALSE) |>
  janitor::clean_names()
if (!inherits(raw$date, "Date")) raw$date <- as.Date(raw$date)
stopifnot(!anyNA(raw$date))

# Use official national totals for national trends, separately from states.
# Exclude the residual "Others" geography from the state-level panel.
selected <- raw |>
  filter(!state %in% c("India Total", "Others")) |>
  filter(
    (category == "Electricity generation" & unit == "GWh" &
       subcategory %in% c("Fuel", "Total")) |
    (category == "Capacity" & unit == "MW" & subcategory == "Fuel")
  ) |>
  select(state, date, category, subcategory, variable, value)

# Combine Jammu and Kashmir with Ladakh to hold the geography consistent.
# Sum only observed values; never interpret an entirely missing group as zero.
panel_long <- selected |>
  mutate(
    state = if_else(
      state %in% c("Jammu and Kashmir", "Ladakh"),
      "Jammu and Kashmir + Ladakh", state
    )
  ) |>
  group_by(state, date, category, subcategory, variable) |>
  summarise(
    value = if (all(is.na(value))) NA_real_ else sum(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    fuel_name = janitor::make_clean_names(variable, allow_dupes = TRUE),
    measure = case_when(
      category == "Capacity" ~ paste0("cap_", fuel_name),
      subcategory == "Total" ~ "gen_total",
      TRUE ~ paste0("gen_", fuel_name)
    )
  ) |>
  select(state, date, measure, value)

# The keys must be unique before pivoting to a wide panel.
stopifnot(!anyDuplicated(panel_long[c("state", "date", "measure")]))

panel <- panel_long |>
  pivot_wider(
    id_cols = c(state, date),
    names_from = measure, values_from = value
  ) |>
  arrange(state, date)

required_columns <- c(
  "gen_total", "gen_coal", "gen_solar", "gen_wind", "cap_coal"
)
stopifnot(all(required_columns %in% names(panel)))

# The coal capacity factor is generation / potential generation.
# GWh -> MWh requires multiplication by 1,000, equivalent to /1,000
# in the denominator when capacity is in MW.
panel <- panel |>
  mutate(
    coal_share = if_else(gen_total > 0, gen_coal / gen_total, NA_real_),
    solar_share = if_else(gen_total > 0, gen_solar / gen_total, NA_real_),
    vre_share = if_else(
      gen_total > 0, (gen_solar + gen_wind) / gen_total, NA_real_
    ),
    coal_cf = if_else(
      !is.na(cap_coal) & cap_coal > 0,
      gen_coal / (cap_coal * lubridate::days_in_month(date) * 24 / 1000),
      NA_real_
    )
  )

# These assertions refer to the Ember release used for this sample.
stopifnot(
  nrow(panel) == 2905,
  n_distinct(panel$state) == 35,
  n_distinct(panel$date) == 83,
  !anyDuplicated(panel[c("state", "date")]),
  all(is.na(panel$coal_share) | between(panel$coal_share, 0, 1)),
  all(is.na(panel$vre_share) | between(panel$vre_share, 0, 1)),
  all(is.na(panel$coal_cf) | between(panel$coal_cf, 0, 1))
)

# Compare total generation to its component fuels on COMPLETE rows only.
# Rows with missing fuels are reported, not silently treated as zero.
fuel_columns <- setdiff(grep("^gen_", names(panel), value = TRUE), "gen_total")
fuel_balance <- panel |>
  mutate(
    n_missing_fuels = rowSums(is.na(pick(all_of(fuel_columns)))),
    fuel_sum = rowSums(pick(all_of(fuel_columns)), na.rm = TRUE),
    balance_gap_gwh = if_else(
      n_missing_fuels == 0, gen_total - fuel_sum, NA_real_
    )
  )

quality_checks <- fuel_balance |>
  summarise(
    state_months = n(),
    incomplete_fuel_rows = sum(n_missing_fuels > 0),
    max_complete_row_gap_gwh = max(abs(balance_gap_gwh), na.rm = TRUE),
    gaps_over_1_gwh = sum(abs(balance_gap_gwh) > 1, na.rm = TRUE)
  )
print(quality_checks)
stopifnot(quality_checks$gaps_over_1_gwh == 0)

missing_report <- panel |>
  summarise(across(all_of(required_columns), ~ round(100 * mean(is.na(.x)), 2))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "missing_pct") |>
  arrange(desc(missing_pct))
print(missing_report)

saveRDS(panel, file.path("data", "processed", "panel.rds"))
cat("\nPrepared", nrow(panel), "state-month observations across",
    n_distinct(panel$state), "regions.\n")

# 3. DESCRIPTIVE ANALYSIS -----------------------------------------------------

# Figure 1: Monthly national generation from coal, solar and wind.
national_gen <- raw |>
  filter(
    state == "India Total", category == "Electricity generation",
    subcategory == "Fuel", unit == "GWh",
    variable %in% c("Coal", "Solar", "Wind")
  ) |>
  select(date, variable, value) |>
  mutate(variable = factor(variable, levels = c("Coal", "Solar", "Wind")))
stopifnot(nrow(national_gen) == 249,
          !anyDuplicated(national_gen[c("date", "variable")]))

fig1 <- ggplot(national_gen, aes(date, value)) +
  geom_line(color = "#2563A6", linewidth = 0.9) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Electricity Generation in India",
    subtitle = "Monthly coal, solar and wind generation, 2019–2025",
    x = NULL, y = "Generation (GWh)",
    caption = "Source: Ember. National grid electricity generation."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
save_figure(fig1, "01_national_generation_trends.png", height = 8)

# Figure 2: National fuel mix, Jan–Nov consistently for all years.
national_total <- raw |>
  filter(
    state == "India Total", category == "Electricity generation",
    subcategory == "Total", unit == "GWh"
  ) |>
  select(date, total_gwh = value)
stopifnot(!anyDuplicated(national_total$date), !anyNA(national_total$total_gwh))

national_mix <- national_gen |>
  mutate(variable = as.character(variable)) |>
  pivot_wider(names_from = variable, values_from = value) |>
  left_join(national_total, by = "date") |>
  filter(lubridate::month(date) <= 11) |>
  group_by(year = lubridate::year(date)) |>
  summarise(
    months = n(),
    coal_gwh = sum(Coal), solar_gwh = sum(Solar),
    wind_gwh = sum(Wind), total_gwh = sum(total_gwh),
    .groups = "drop"
  ) |>
  mutate(
    coal_share = 100 * coal_gwh / total_gwh,
    solar_share = 100 * solar_gwh / total_gwh,
    wind_share = 100 * wind_gwh / total_gwh,
    vre_share = solar_share + wind_share,
    other_share = 100 - coal_share - solar_share - wind_share
  )
stopifnot(all(national_mix$months == 11),
          all(between(national_mix$other_share, 0, 100)))

mix_long <- national_mix |>
  select(year, coal_share, solar_share, wind_share, other_share) |>
  pivot_longer(-year, names_to = "source", values_to = "share") |>
  mutate(source = factor(
    source,
    levels = c("coal_share", "solar_share", "wind_share", "other_share"),
    labels = c("Coal", "Solar", "Wind", "Other sources")
  ))

fig2 <- ggplot(mix_long, aes(factor(year), share, fill = source)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Coal" = "#4B5563", "Solar" = "#F59E0B",
    "Wind" = "#16A6A0", "Other sources" = "#CBD5E1"
  )) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = seq(0, 100, 20), limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "India's Electricity Generation Mix",
    subtitle = "Fuel shares, January–November of each year",
    x = NULL, y = "Share of electricity generation", fill = NULL,
    caption = "Source: Ember. Grid generation excluding captive power. Other sources include hydro, gas, nuclear and bioenergy."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top", panel.grid.major.x = element_blank())
save_figure(fig2, "02_generation_mix.png")
readr::write_csv(national_mix, file.path("tables", "national_generation_mix.csv"))

# Figure 3: State-level solar + wind share (same partial-year window).
regional_mix <- panel |>
  mutate(year = lubridate::year(date), month = lubridate::month(date)) |>
  filter(month <= 11) |>
  group_by(state, year) |>
  summarise(
    months = n(),
    solar_gwh = sum(gen_solar), wind_gwh = sum(gen_wind),
    total_gwh = sum(gen_total), .groups = "drop"
  ) |>
  mutate(vre_share_pct = if_else(
    total_gwh > 0,
    100 * (solar_gwh + wind_gwh) / total_gwh,
    NA_real_
  ))
stopifnot(all(regional_mix$months == 11))

state_order <- regional_mix |>
  filter(year == 2025) |>
  arrange(vre_share_pct) |>
  pull(state)

fig3 <- regional_mix |>
  mutate(state = factor(state, levels = state_order)) |>
  ggplot(aes(factor(year), state, fill = vre_share_pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "#F1F5F9", high = "#0F766E", na.value = "grey85",
    limits = c(0, 100), name = "Solar + wind\nshare (%)"
  ) +
  labs(
    title = "Solar and Wind Shares Across Indian Regions",
    subtitle = "January–November of each year",
    x = NULL, y = NULL,
    caption = "Source: Ember. Shares can be high where local total generation is small."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank(), axis.text.y = element_text(size = 8))
save_figure(fig3, "03_regional_renewables_heatmap.png", height = 11)

# Figure 4: Absolute VRE production (TWh) rather than percentages.
top_vre_states <- regional_mix |>
  filter(year == 2025) |>
  mutate(vre_twh = (solar_gwh + wind_gwh) / 1000) |>
  slice_max(vre_twh, n = 10, with_ties = FALSE)

fig4 <- ggplot(
  top_vre_states,
  aes(vre_twh, forcats::fct_reorder(state, vre_twh))
) +
  geom_col(fill = "#0F766E", width = 0.7) +
  geom_text(aes(label = round(vre_twh, 1)), hjust = -0.2, size = 3.5) +
  scale_x_continuous(
    labels = scales::comma, expand = expansion(mult = c(0, 0.13))
  ) +
  labs(
    title = "Solar and Wind Generation Across Indian Regions",
    subtitle = "Top 10 regions, January–November 2025",
    x = "Generation (TWh)", y = NULL,
    caption = "Source: Ember. State-level grid electricity generation."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
save_figure(fig4, "04_renewable_generation_by_state.png")

# Figure 5: Monthly seasonality, using only complete calendar years.
seasonal_data <- national_gen |>
  mutate(year = lubridate::year(date), month = lubridate::month(date)) |>
  filter(year >= 2019, year <= 2024) |>
  group_by(variable, year) |>
  mutate(monthly_index = 100 * value / mean(value)) |>
  ungroup() |>
  group_by(variable, month) |>
  summarise(
    seasonal_index = mean(monthly_index), years = n(), .groups = "drop"
  )
stopifnot(nrow(seasonal_data) == 36, all(seasonal_data$years == 6))

fig5 <- ggplot(seasonal_data, aes(month, seasonal_index)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey60") +
  geom_line(color = "#2563A6", linewidth = 0.9) +
  geom_point(color = "#2563A6", size = 2) +
  facet_wrap(~ variable, ncol = 1) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Seasonality of Electricity Generation in India",
    subtitle = "Mean monthly index across complete years, 2019–2024",
    x = NULL, y = "Index (annual monthly mean = 100)",
    caption = "Source: Ember; author's calculations."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
save_figure(fig5, "05_generation_seasonality.png", height = 8)

# 4. REGRESSION DESIGN --------------------------------------------------------

# Dependent variable: coal capacity utilization (%) in a state-month.
# Main explanatory variable: solar + wind share of generation (%).
# State and monthly time fixed effects; standard errors clustered by state.
reg_sample <- panel |>
  mutate(
    coal_cf_pp = 100 * coal_cf,
    vre_share_pp = 100 * vre_share
  ) |>
  filter(!is.na(coal_cf_pp), !is.na(vre_share_pp), cap_coal > 0)

# Delhi and Lakshadweep have sparse, zero-coal-generation eligible rows.
# The extended-sample model below retains them as a sensitivity check.
core_sample <- reg_sample |>
  filter(!state %in% c("Delhi", "Lakshadweep")) |>
  mutate(
    solar_share_pp = 100 * gen_solar / gen_total,
    wind_share_pp = 100 * gen_wind / gen_total,
    solar_twh = gen_solar / 1000,
    wind_twh = gen_wind / 1000,
    time_index = 12 * (lubridate::year(date) - 2019) + lubridate::month(date)
  )
stopifnot(
  nrow(core_sample) == 1377,
  n_distinct(core_sample$state) == 17,
  n_distinct(core_sample$date) == 81,
  !anyNA(core_sample$solar_share_pp),
  !anyNA(core_sample$wind_share_pp),
  all(abs(core_sample$solar_share_pp + core_sample$wind_share_pp -
            core_sample$vre_share_pp) < 1e-8)
)
cat("\nRegression sample:", nrow(core_sample), "observations in",
    n_distinct(core_sample$state), "states.\n")

model_fe <- fixest::feols(
  coal_cf_pp ~ vre_share_pp | state + date,
  data = core_sample, vcov = ~ state
)
model_all <- fixest::feols(
  coal_cf_pp ~ vre_share_pp | state + date,
  data = reg_sample, vcov = ~ state
)
model_split <- fixest::feols(
  coal_cf_pp ~ solar_share_pp + wind_share_pp | state + date,
  data = core_sample, vcov = ~ state
)
# In this reparameterization, the wind coefficient is beta_wind-beta_solar.
model_difference <- fixest::feols(
  coal_cf_pp ~ vre_share_pp + wind_share_pp | state + date,
  data = core_sample, vcov = ~ state
)
# Absolute generation levels avoid putting national/regional total generation
# in the explanatory variables' denominator.
model_levels <- fixest::feols(
  coal_cf_pp ~ solar_twh + wind_twh | state + date,
  data = core_sample, vcov = ~ state
)
model_trend_share <- fixest::feols(
  coal_cf_pp ~ vre_share_pp | state[time_index] + date,
  data = core_sample, vcov = ~ state
)
model_trend_levels <- fixest::feols(
  coal_cf_pp ~ solar_twh + wind_twh | state[time_index] + date,
  data = core_sample, vcov = ~ state
)

cat("\nMain and alternative regression specifications:\n")
print(fixest::etable(
  model_fe, model_all, model_split, model_levels,
  model_trend_share, model_trend_levels, digits = 3
))
cat("\nTest of the difference between solar and wind coefficients:\n")
print(fixest::coeftable(model_difference))

# 5. FINITE-CLUSTER INFERENCE: WILD CLUSTER BOOTSTRAP -------------------------

# With 17 state clusters, also assess inference via a wild cluster bootstrap.
# Factor-based FE avoid a Date-processing error observed in an earlier run.
boot_data <- core_sample |>
  mutate(state_boot = factor(state), date_boot = factor(date))

model_boot <- fixest::feols(
  coal_cf_pp ~ vre_share_pp | state_boot + date_boot,
  data = boot_data, vcov = ~ state_boot
)
# An algebraically equivalent trend specification for the bootstrap package.
model_trend_boot <- fixest::feols(
  coal_cf_pp ~ vre_share_pp +
    i(state_boot, time_index, ref = "Andhra Pradesh") |
    state_boot + date_boot,
  data = boot_data, vcov = ~ state_boot
)

stopifnot(
  isTRUE(all.equal(
    unname(coef(model_fe)["vre_share_pp"]),
    unname(coef(model_boot)["vre_share_pp"]), tolerance = 1e-8
  )),
  isTRUE(all.equal(
    unname(coef(model_trend_share)["vre_share_pp"]),
    unname(coef(model_trend_boot)["vre_share_pp"]), tolerance = 1e-6
  ))
)

# 9,999 Rademacher bootstrap replications, clustered by state.
set.seed(2026)
dqrng::dqset.seed(2026)
boot_final <- fwildclusterboot::boottest(
  model_boot, param = "vre_share_pp", clustid = "state_boot",
  B = 9999, sampling = "dqrng", conf_int = FALSE
)
set.seed(2026)
dqrng::dqset.seed(2026)
boot_trends_final <- fwildclusterboot::boottest(
  model_trend_boot, param = "vre_share_pp", clustid = "state_boot",
  B = 9999, sampling = "dqrng", conf_int = FALSE
)
cat("\nBaseline bootstrap:\n")
print(summary(boot_final))
cat("\nState-trend bootstrap:\n")
print(summary(boot_trends_final))

# 6. LEAVE-ONE-STATE-OUT ROBUSTNESS ------------------------------------------

loo_results <- purrr::map_dfr(
  sort(unique(core_sample$state)),
  function(excluded_state) {
    sample_reduced <- core_sample |>
      filter(state != excluded_state)
    model_loo <- fixest::feols(
      coal_cf_pp ~ vre_share_pp | state + date,
      data = sample_reduced, vcov = ~ state
    )
    ct <- fixest::coeftable(model_loo)
    tibble(
      excluded_state = excluded_state,
      beta = unname(ct["vre_share_pp", "Estimate"]),
      std_error = unname(ct["vre_share_pp", "Std. Error"]),
      p_value = unname(ct["vre_share_pp", "Pr(>|t|)"])
    )
  }
)

loo_summary <- loo_results |>
  summarise(
    n_models = n(), min_beta = min(beta), max_beta = max(beta),
    all_negative = all(beta < 0),
    n_clustered_p_below_0_05 = sum(p_value < 0.05)
  )
print(loo_summary)

fig6 <- loo_results |>
  mutate(excluded_state = forcats::fct_reorder(excluded_state, beta)) |>
  ggplot(aes(beta, excluded_state)) +
  geom_vline(
    xintercept = unname(coef(model_fe)["vre_share_pp"]),
    linetype = "dashed", color = "#D97757", linewidth = 0.8
  ) +
  geom_point(color = "#2563A6", size = 2.8) +
  labs(
    title = "Sensitivity to Individual States",
    subtitle = "Leave-one-state-out estimates of the VRE coefficient",
    x = "Estimated coefficient", y = "Excluded state",
    caption = "Dashed line: full-sample estimate. Source: Ember; author's calculations."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())
save_figure(fig6, "06_leave_one_state_out.png")

# 7. EXPORT THE FINAL REGRESSION RESULTS -------------------------------------

regression_summary <- tibble(
  model = c("Two-way fixed effects", "Two-way FE + state trends"),
  coefficient = c(
    unname(coef(model_fe)["vre_share_pp"]),
    unname(coef(model_trend_share)["vre_share_pp"])
  ),
  clustered_se = c(
    unname(fixest::se(model_fe)["vre_share_pp"]),
    unname(fixest::se(model_trend_share)["vre_share_pp"])
  ),
  bootstrap_p = c(
    fwildclusterboot::pval(boot_final),
    fwildclusterboot::pval(boot_trends_final)
  ),
  observations = c(nobs(model_fe), nobs(model_trend_share)),
  r_squared = c(
    as.numeric(fixest::fitstat(model_fe, "r2")),
    as.numeric(fixest::fitstat(model_trend_share, "r2"))
  )
)
print(regression_summary)

readr::write_csv(
  regression_summary, file.path("tables", "main_regression_results.csv")
)
readr::write_csv(
  loo_results, file.path("tables", "leave_one_state_out.csv")
)
saveRDS(
  list(
    baseline = model_fe, extended = model_all,
    separate_shares = model_split, difference = model_difference,
    generation_levels = model_levels, state_trends = model_trend_share,
    state_trend_levels = model_trend_levels
  ),
  file.path("data", "processed", "regression_models.rds")
)

stopifnot(
  file.exists(file.path("tables", "main_regression_results.csv")),
  file.exists(file.path("figures", "06_leave_one_state_out.png"))
)
cat("\nCompleted. Figures, tables and processed data saved successfully.\n")

# Limitations:
# - Observational FE associations; no exogenous identification strategy.
# - VRE shares and coal utilization may respond jointly to demand and dispatch.
# - VRE shares are compositional, making level-based checks important.
# - Coal-capacity records are incomplete for some months / regions.
# - Only 17 baseline state clusters; wild-bootstrap p-values are reported.
# - Absolute-generation specifications provide less consistent evidence.
