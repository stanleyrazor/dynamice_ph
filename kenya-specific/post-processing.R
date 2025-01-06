

pacman::p_load(dplyr, ggplot2, patchwork, tidyr, stringr, purrr, gt)

# data
d1 <- read.csv("central_burden_estimate/central_burden_estimate_nomcv.csv") |>
  select(country, year, age,
         starts_with("cases"),
         starts_with("death"), dalys) |>
  mutate(cases = cases0d + cases1d + cases2d,
         deaths = deaths0d + deaths1d + deaths2d,
         strategy = "No MCV") |>
  select(strategy, country, year, age, cases, deaths, dalys)

d2 <- read.csv("central_burden_estimate/central_burden_estimate_mcv1.csv") |>
  select(country, year, age,
         starts_with("cases"),
         starts_with("death"), dalys) |>
  mutate(cases = cases0d + cases1d + cases2d,
         deaths = deaths0d + deaths1d + deaths2d,
         strategy = "MCV1") |>
  select(strategy, country, year, age, cases, deaths, dalys)

d3 <- read.csv("central_burden_estimate/central_burden_estimate_mcv1-mcv2.csv") |>
  select(country, year, age,
         starts_with("cases"),
         starts_with("death"), dalys) |>
  mutate(cases = cases0d + cases1d + cases2d,
         deaths = deaths0d + deaths1d + deaths2d,
         strategy = "MCV1 & MCV2") |>
  select(strategy, country, year, age, cases, deaths, dalys)

# f1 <- read.csv("~/Documents/GitHub/dynamice/central_burden_estimate/county/Portnoy/cbe_no-vacc_syn_linearVE_Portnoy.csv")
# f2 <- f1 |>
#   select(country, year, age, cases, deaths, dalys)
#
# list(f3 = f2, d3 = d2) |>
#   lapply(\(data) {
#     data |>
#       group_by(country, year) |>
#       reframe(cases = sum(cases),
#               deaths = sum(deaths),
#               dalys = sum(dalys))
#   }) |>
#   list2env(envir = .GlobalEnv)


# -------------------------------------------------------------------------

d4 <- bind_rows(d1, d2, d3)

## 2013 to 2023
d5 <- d4 |>
  mutate(period = ifelse(year %in% 2013:2023, "Actual", "Projection")) |>
  group_by(country, period, strategy) |>
  reframe(cases = sum(cases),
          deaths = sum(deaths),
          dalys = sum(dalys))

averted_cases <- d5 |>
  select(country, period, strategy, cases) |>
  pivot_wider(names_from = strategy, values_from = cases) |>

  mutate(
    `Averted due to MCV1` = -1 * (MCV1 - `No MCV`),
    `Averted due to MCV2` = -1 * (`MCV1 & MCV2` -MCV1)
  ) |>
  select(county = country, Period = period, starts_with("Averted"))

averted_deaths <- d5 |>
  select(country, period, strategy, deaths) |>
  pivot_wider(names_from = strategy, values_from = deaths) |>

  mutate(
    `Averted due to MCV1` = -1 * (MCV1 - `No MCV`),
    `Averted due to MCV2` = -1 * (`MCV1 & MCV2` -MCV1)
  ) |>
  select(county = country, Period = period, starts_with("Averted"))

averted_dalys <- d5 |>
  select(country, period, strategy, dalys) |>
  pivot_wider(names_from = strategy, values_from = dalys) |>

  mutate(
    `Averted due to MCV1` = -1 * (MCV1 - `No MCV`),
    `Averted due to MCV2` = -1 * (`MCV1 & MCV2` -MCV1)
  ) |>
  select(county = country, Period = period, starts_with("Averted"))


# National tables ---------------------------------------------------------

national_actual_burden <- d5 |>
  group_by(strategy, period) |>
  reframe(across(where(is.numeric), sum)) |>
  mutate(across(where(is.numeric), \(x) round(x, 0)))

national_averted_burden <- bind_rows(averted_cases |> mutate(burden = "Cases")
          , averted_deaths |> mutate(burden = "Deaths")
          , averted_dalys |> mutate(burden = "DALYs")
          ) |>
  group_by(burden, Period) |>
  reframe(across(where(is.numeric), sum)) |>
  mutate(across(where(is.numeric), \(x) round(x, 0)))


# -------------------------------------------------------------------------

# national actual burden
national_actual_burden |>
  pivot_wider(
    names_from = period,
    values_from = c(cases, deaths, dalys)
  ) |>
  mutate(strategy = factor(strategy, levels = c("No MCV", "MCV1", "MCV1 & MCV2"))) |>
  arrange(strategy) |>
  gt(rowname_col = "strategy") |>
  tab_spanner(
    label = "Actual (2013-2023)",
    columns = c(cases_Actual, deaths_Actual, dalys_Actual)
  ) |>
  tab_spanner(
    label = "Projection (2024-2030)",
    columns = c(cases_Projection, deaths_Projection, dalys_Projection)
  ) |>
  cols_label(
    cases_Actual = "Cases",
    deaths_Actual = "Deaths",
    dalys_Actual = "DALYs",
    cases_Projection = "Cases",
    deaths_Projection = "Deaths",
    dalys_Projection = "DALYs"
  ) |>
  fmt_number(
    columns = where(is.numeric),
    use_seps = TRUE,
    decimals = 0
  ) |>
  tab_options(
    table.font.names = "Times New Roman"
  )


# national averted burden
national_averted_burden |>
  pivot_wider(
    names_from = Period,
    values_from = c(`Averted due to MCV1`, `Averted due to MCV2`)
  ) |>
  gt(rowname_col = "burden") |>
  tab_spanner(
    label = "Actual (2013-2023)",
    columns = c(`Averted due to MCV1_Actual`, `Averted due to MCV2_Actual`)
  ) |>
  tab_spanner(
    label = "Projection (2024-2030)",
    columns = c(`Averted due to MCV1_Projection`, `Averted due to MCV2_Projection`)
  ) |>
  cols_label(
    `Averted due to MCV1_Actual` = "Averted due to MCV1",
    `Averted due to MCV2_Actual` = "Averted due to MCV2",
    `Averted due to MCV1_Projection` = "Averted due to MCV1",
    `Averted due to MCV2_Projection` = "Averted due to MCV2"
  ) |>
  fmt_number(
    columns = where(is.numeric),
    use_seps = TRUE,
    decimals = 0
  ) |>
  tab_options(
  table.font.names = "Times New Roman"
)

national_ts <- d4 |>
  group_by(strategy, year) |>
  reframe(across(c(cases, deaths, dalys), sum)) |>
  setNames(c("strategy", "year", "Cases", "Deaths", "DALYs")) |>
  pivot_longer(-c(strategy, year)) |>
  ggplot() +
  geom_line(aes(x = year, y = value, col = strategy)) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  labs(title = "National Measles Health Burden",
       x = NULL, y = "Health burden",
       col = "Vaccine:") +
  theme_bw(base_line_size = 0,
           base_family = "Times New Roman") +
  theme(legend.position = "bottom")

ggsave("img/national_timeseries.png", plot = national_ts,
       width = 8, height = 8, units = "in", dpi = 300)


# County graphs' ----------------------------------------------------------

c1 <- d5 |>
  pivot_longer(c(cases, deaths, dalys)) |>
  mutate(strategy = factor(strategy, levels = c("No MCV", "MCV1", "MCV1 & MCV2"))) |>
  filter(name == "cases") |>
  ggplot(aes(x = reorder(country, value), y = value)) +
  geom_col() +
  facet_wrap( ~ strategy) +
  labs(title = "Measles Health Burden: Cases",
       x = NULL, y = "Health burden") +
  coord_flip() +
  theme_bw(base_line_size = 0, base_family = "Times New Roman")

c2 <- d5 |>
  pivot_longer(c(cases, deaths, dalys)) |>
  mutate(strategy = factor(strategy, levels = c("No MCV", "MCV1", "MCV1 & MCV2"))) |>
  filter(name == "deaths") |>
  ggplot(aes(x = reorder(country, value), y = value)) +
  geom_col() +
  facet_wrap( ~ strategy) +
  labs(title = "Measles Health Burden: Deaths",
       x = NULL, y = "Health burden") +
  coord_flip() +
  theme_bw(base_line_size = 0, base_family = "Times New Roman")

c3 <- d5 |>
  pivot_longer(c(cases, deaths, dalys)) |>
  mutate(strategy = factor(strategy, levels = c("No MCV", "MCV1", "MCV1 & MCV2"))) |>
  filter(name == "dalys") |>
  ggplot(aes(x = reorder(country, value), y = value)) +
  geom_col() +
  facet_wrap( ~ strategy) +
  labs(title = "Measles Health Burden: DALYs",
       x = NULL, y = "Health burden") +
  coord_flip() +
  theme_bw(base_line_size = 0, base_family = "Times New Roman")

ggsave("img/county_healthburden_cases.png", plot = c1,
width = 14, height = 8, units = "in", dpi = 300)

ggsave("img/county_healthburden_deaths.png", plot = c2,
       width = 14, height = 8, units = "in", dpi = 300)

ggsave("img/county_healthburden_dalys.png", plot = c3,
       width = 14, height = 8, units = "in", dpi = 300)

ts_counties <- d4 |>
  group_by(country, year, strategy) |>
  reframe(cases = sum(cases),
          deaths = sum(deaths),
          dalys = sum(dalys)) |>
  ggplot(aes(x = year, y = cases)) +
  geom_line(aes(col = strategy)) +
  labs(title = "Health burden over time",
       x = NULL, y = "Health Burden: Cases",
       col = "Vaccine") +
  theme_bw(base_line_size = 0, base_family = "Times New Roman") +
  theme(legend.position = "bottom") +
  facet_wrap(~country, scales = "free")

ggsave("img/county_timeseries.png", plot = ts_counties,
       width = 14, height = 8, units = "in", dpi = 300)

# -------------------------------------------------------------------------

# Transform data into a wide format
d5_wide <- d5 %>%
  pivot_wider(
    names_from = period,
    values_from = c(cases, deaths, dalys)
  )

# Create the gt table
d5_gt <- d5_wide %>%
  gt(rowname_col = "country", groupname_col = "strategy") %>%
  tab_spanner(
    label = "Actual",
    columns = c(cases_Actual, deaths_Actual, dalys_Actual)
  ) %>%
  tab_spanner(
    label = "Projection",
    columns = c(cases_Projection, deaths_Projection, dalys_Projection)
  ) %>%
  cols_label(
    cases_Actual = "Cases",
    deaths_Actual = "Deaths",
    dalys_Actual = "DALYs",
    cases_Projection = "Cases",
    deaths_Projection = "Deaths",
    dalys_Projection = "DALYs"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    use_seps = TRUE,
    decimals = 0
  ) %>%
  tab_options(
    table.font.names = "Times New Roman"
  )

# View the table
d5_gt
gt::gtsave(d5_gt, filename = "img/cumulative_health_burden.docx")



# Age distribution --------------------------------------------------------

agedist <- d4 |>
  mutate(group = case_when(age >= 3 ~ "3+",
                           TRUE ~ as.character(age))) |>
  group_by(strategy, year, group) |>
  reframe(across(c(cases, deaths, dalys), sum)) |>
  setNames(c("strategy", "year", "group", "Cases", "Deaths", "DALYs")) |>
  pivot_longer(-c(strategy, year, group)) |>
  group_by(strategy, year, name) |>
  mutate(value = value / sum(value)) |>
  ungroup() |>
  filter(strategy == "MCV1 & MCV2") |>
  ggplot(aes(x = year, y = value, fill = group)) +
  geom_area() +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(. ~ name, nrow = 2) +
  labs(
    title = "Age distribution of cases",
    x = "Year",
    y = "Percentage distribution",
    fill = NULL
  ) +
  scale_fill_manual(values = c("#27aae1", "#fe7501", "#084887", "#efca08")) +
  theme_minimal(base_line_size = 0, base_family = "Times New Roman") +
  theme(legend.position = "bottom")


ggsave("img/age_dist.png", plot = agedist,
       width = 8, height = 8, units = "in", dpi = 300)


