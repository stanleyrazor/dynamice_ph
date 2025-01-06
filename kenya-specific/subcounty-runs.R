# run_scenario
# execute simulations for different vaccination scenarios
# update: 2023/04/07

rm(list = ls())

# load libraries, data, functions
pacman::p_load (data.table, dplyr, stringr, Rcpp, doParallel)

# rcpp codes
sourceCpp("rcpp/rcpp_spinup.cpp")
sourceCpp("rcpp/rcpp_vaccine_oney.cpp")

load(file = "data/data_pop.rda")
load(file = "data/data_cfr_portnoy_21.rda")
load(file = "data/data_contact_syn.rda")
load(file = "data/data_r0.rda")
load(file = "data/data_timeliness.rda")
load(file = "data/data_lexp_remain.rda")
load(file = "data/data_template.rda")

# list countries for analysis
county_codes <- read.csv("kenya-specific/data/county_codes.csv") |>
  mutate(name = str_to_title(name))
use_county_code <- county_codes |> pull(code)
use_county_name <- county_codes |> pull(name)

eva_code <- eva_countries <- use_county_code[-1]
eva_name <- county_codes |> filter(code %in% eva_code) |> pull(name)

# -------------------------------------------------------------------------

# MODIFYING THE CONTACT MATRIX TO ONLY BE KENYAN
data_contact_syn <- list(data_contact_syn[["KEN"]])
names(data_contact_syn) = "KEN"

# adding county level data to the timeliness dataset
county_timeliness <- expand.grid(
  country_code = use_county_code[-1],
  age = 0:152
) |>
  mutate(timeliness = ifelse(age < 52*(9/12), 0, 1)) |>
  bind_rows(
    expand.grid(country_code = use_county_code[-1],
                prop_final_cov = 0.83908) # this is kenya'sprop_final_coverage: not really sure what it is, and whether it is even used.
  )

data_timeliness <- bind_rows(data_timeliness, county_timeliness)

## updating the data template to have kenyan counties
data_template <- data_template |>
  bind_rows(
    expand.grid(
      disease = "Measles",
      year = 2000:2100,
      age = 0:99,
      country = use_county_code[-1],
      cohort_size = NA,
      deaths = NA,
      cases = NA,
      dalys = NA
    ) |>
      merge(data.frame(country = use_county_code[-1],
                       country_name = use_county_name[-1]),
            by = "country")
  )

## updating the population data
kenya_pop <- data_pop |> filter(country == "Kenya")

# population data
county_pop <- (rKenyaCensus::V3_T2.3) |>
  filter(Age %in% c(0:99, "100+") & SubCounty == "ALL") |>
  select(county = County, age = Age, pop = Total) |>
  mutate(age = ifelse(age == "100+", 100, age)) |>
  group_by(age) |>
  mutate(prop = pop / sum(pop),
         county = str_to_lower(county),
         county = ifelse(county == "taita/ taveta", "taita taveta", county),
         county = str_to_title(county)) |>
  ungroup() |>
  select(county, age, prop)

county_pop <- expand.grid(county = county_pop |> pull(county) |> unique(),
                          year = 2000:2030) |>
  merge(kenya_pop |> select(year, age_from, age_to, value), by = "year") |>
  merge(county_pop |> select(county, age_from = age, prop), by = c("county", "age_from")) |>
  mutate(value = value * prop,
         gender = "both") |>
  merge(county_codes |> select(county = name, country_code = code), by = "county") |>
  select(country_code, country = county, age_from, age_to, year, gender, value)

data_pop <- bind_rows(county_pop,
                      kenya_pop |> select(colnames(county_pop)))

## Modifying the life expectancy data
temp_lexp_data <- NULL
for (i in 1:length(eva_countries)) {
  temp_lexp_data <- bind_rows(temp_lexp_data,
                              data_lexp_remain |> filter(country == "Kenya") |>
                                mutate(country_code = eva_code[i],
                                       country = eva_name[i]))
}
data_lexp_remain <- temp_lexp_data

## modifying the CFR Portnoy dataset
temp_cfr_data <- NULL
for (i in 1:length(eva_countries)) {
  temp_cfr_data <- bind_rows(temp_cfr_data,
                             data_cfr_portnoy_21 |> filter(country == "KEN") |>
                               mutate(country_code = eva_code[i],
                                      country = eva_name[i]))
}
data_cfr_portnoy_21 <- temp_cfr_data


# -------------------------------------------------------------------------

# # check percentage of birth cohorts
# sum(data_pop [age_from == 0 & year %in% 2010:2019 & country_code %in% eva_ctries, value])/
#   sum(data_pop [age_from == 0 & year %in% 2010:2019, value])
# # 52.7%

# add data for country-specific age at vaccination
# https://immunizationdata.who.int/pages/schedule-by-disease/measles.html
# input monthly age and then covert to weekly age for dynaMICE structure
data_vage <- data.table(country_code = eva_countries,
                        mcv1 = (9 / 12) * 52,
                        mcv2 = (18 / 12) * 52)

# adjust to yearly age for those >= 3 years old
data_vage [mcv2 > 52*3, mcv2 := 52*3 + floor(mcv2/52-2)]

# update timeliness data for countries not giving MCV1 to 39 week old (9 months)
# MUTED because kenya gives it at that time
# for (ictry in data_vage [mcv1 != 39, country_code]){
#   c_age <- data_vage [country_code == ictry, mcv1]
#   data_timeliness [country_code == ictry & !is.na(age), timeliness := ifelse (age < c_age, 0, 1)]
# }

# assume a fixed R0 for the central run
# median R0 for least developed countries, vaccine era, from Guerra et al. (2017)
adj.fixR0 = 14.25156  # NA
# {if (!is.na(adj.fixR0)){
#   data_r0 <- copy (data_r0 [country_code %in% eva_countries]) [ , r0 := adj.fixR0]
# }}

data_r0 <- # data_r0 |>
  bind_rows(data.frame(country = use_county_name,
                       country_code = use_county_code,
                       r0 = adj.fixR0))

# load functions
source ("R/logs.R")
source ("R/utils.R")
source ("R/functions_rcpp.R")

# set up variables
var <- list (
  # vaccine coverage
  vaccine_coverage_folder           = "coverage/",
  coverage_prefix                   = "coverage",
  touchstone                        = "_",
  vaccine_coverage_subfolder        = "scenarios/kenya-county/",

  # disease burden
  burden_estimate_folder            = "central_burden_estimate/",

  # log file name
  log_name                          = "test_log",

  # country iso3 codes
  countries                         = eva_countries
)

# prepare coverage inputs
vac_strategies <- c("nomcv",                 # (1) no vaccination
                    "mcv1",                  # (2) MCV1 only
                    "mcv1-mcv2",             # (3) MCV1 + MCV2
                    "mcv1-mcv2-sia",         # (4) MCV1 + MCV2 + SIA
                    "mcv1-sia",              # (5) MCV1 + SIA
                    "mcv1-mcv2alt1",         # (6) MCV1 + MCV2(early intro, fast rollout)
                    "mcv1-mcv2alt1-sia",     # (7) MCV1 + MCV2(early intro, fast rollout) + SIA
                    "mcv1-mcv2-siaalt1",     # (8) MCV1 + MCV2 + SIA(zero dose first)
                    "mcv1-mcv2-siaalt2",     # (9) MCV1 + MCV2 + SIA(already vaccinated first)
                    "mcv1-siaalt1",          # (10) MCV1 + SIA(zero dose first)
                    "mcv1-siaalt2",          # (11) MCV1 + SIA(already vaccinated first)
                    "mcv1-mcv2alt1-siaalt1", # (12) MCV1 + MCV2(early intro) + SIA(zero dose first)
                    "mcv1-mcv2alt2"          # (13) MCV1 + MCV2(early intro, gradual rollout)
)

# set SIAs implementation method for each scenario
# additional modifications needed in Rcpp for extent-specific approaches
# ex. random reach for subnational campaign & 7.7% at national level for national and rollover-nat campaigns
# 0: no SIA
# 1: random reach (baseline assumption)
# 2: 7.7% less likely to be reached at national level
# 3: zero-dose first
# 4: already-vaccinated first
# 5: 7.7% less likely to be reached at subnational level
set_sia         <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0)

# set routine vaccination parameters to distinguish between scenarios
# 0: no routine MCV
# 1: MCV1 only
# 2: MCV1 + MCV2
set_vaccination <- c (0, 1, 2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2)

# prepare coverage input data - update when the data are changed
adj.covfiles <- FALSE
{if(adj.covfiles){
  # generate 2 coverage input files for routine and SIA vaccination
  for (index in 1:length(vaccine_strategies)){
    create_vaccine_coverage_routine_sia (
      vaccine_coverage_folder    = var$vaccine_coverage_folder,
      vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
      coverage_prefix            = var$coverage_prefix,
      touchstone                 = var$touchstone,
      scenario_name              = vac_strategies [index]
    )
  }
}}

# create a folder for burden estimates results
dir.create (file.path (paste0 (getwd(), "/", var$burden_estimate_folder)), recursive = T)

# removing files
rm(temp_cfr_data, temp_lexp_data, kenya_pop, county_timeliness, county_pop)

# -------------------------------------------------------------------------


# analyse burden estimates under different R0
set_sia <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0) # set back to the oringinal assumptions
ir0 <- 14.25156

# vary R0 values
data_r0 <- copy (setDT(data_r0)) [ , r0 := ir0]
source ("R/functions_rcpp.R")
sel_scns <- c(1,2,3)

for (index in sel_scns) {

  scenario_name   <- vac_strategies [index]
  scenario_number <- sprintf ("scenario%02d", index)

  message(scenario_name)

  # run model and estimate cases
  burden_estimate_file <- runScenario_rcpp (
    vaccine_coverage_folder    = var$vaccine_coverage_folder,
    vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
    coverage_prefix            = var$coverage_prefix,
    scenario_name              = scenario_name,
    save_scenario              = scenario_number,
    burden_estimate_folder     = var$burden_estimate_folder,
    log_name                   = var$log_name,
    countries                  = var$countries,
    vaccination                = set_vaccination [index],
    using_sia                  = set_sia [index],
    sim_years                  = 2013:2030
  )

  # separately estimate deaths
  burden_estimate_file <- paste0 ("central_burden_estimate_",
                                  scenario_name, ".csv")

  # merge outputs into csv files
  get_burden_estimate (vaccine_coverage_folder    = var$vaccine_coverage_folder,
                       vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
                       scenario_name              = vac_strategies [index],
                       save_scenario              = sprintf ("scenario%02d", index),
                       burden_estimate_folder     = var$burden_estimate_folder,
                       log_name                   = var$log_name,
                       vaccination                = set_vaccination [index],
                       using_sia                  = set_sia         [index],
                       folder_date                = Sys.Date() |> str_replace_all("-", ""),  # select the correct folder
                       sim_years                  = 2013:2030)

  # rename file
  # file.rename (from = paste0 (var$burden_estimate_folder, burden_estimate_file),
  #              to = paste0 (var$burden_estimate_folder,
  #                           "central_burden_estimate_",
  #                           scenario_name,
  #                           "_r0-", ir0, ".csv"))
}

# -------------------------------------------------------------------------

cl <- makeCluster(3)
registerDoParallel(cl)

foreach(i = 1:length(sel_scns),
        .packages = c("data.table", "dplyr",
                      "dplyr", "stringr", "Rcpp"),
        .combine = c,
        .noexport = c("rcpp_spinup", "rcpp_vaccine_oney"),
        .verbose = T) %dopar% {

          sourceCpp("rcpp/rcpp_spinup.cpp")
          sourceCpp("rcpp/rcpp_vaccine_oney.cpp")

          index <- sel_scns[i]
          scenario_name   <- vac_strategies [index]
          scenario_number <- sprintf ("scenario%02d", index)

          message(scenario_name)

          # run model and estimate cases
          burden_estimate_file <- runScenario_rcpp (
            vaccine_coverage_folder    = var$vaccine_coverage_folder,
            vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
            coverage_prefix            = var$coverage_prefix,
            scenario_name              = scenario_name,
            save_scenario              = scenario_number,
            burden_estimate_folder     = var$burden_estimate_folder,
            log_name                   = var$log_name,
            countries                  = var$countries,
            vaccination                = set_vaccination [index],
            using_sia                  = set_sia [index],
            sim_years                  = 2013:2030
          )

          # separately estimate deaths
          burden_estimate_file <- paste0 ("central_burden_estimate_",
                                          scenario_name, ".csv")

          # merge outputs into csv files
          get_burden_estimate (vaccine_coverage_folder    = var$vaccine_coverage_folder,
                               vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
                               scenario_name              = vac_strategies [index],
                               save_scenario              = sprintf ("scenario%02d", index),
                               burden_estimate_folder     = var$burden_estimate_folder,
                               log_name                   = var$log_name,
                               vaccination                = set_vaccination [index],
                               using_sia                  = set_sia         [index],
                               folder_date                = Sys.Date() |> str_replace_all("-", ""),  # select the correct folder
                               sim_years                  = 2013:2030)

          burden_estimate_file;
        }

beepr::beep(2)
stopCluster(cl)

# -------------------------------------------------------------------------

#
#
# # move files to a specified folder
# res_files <- list.files (var$burden_estimate_folder)
# dir.create (paste0 ("previous_res/20230401/siareach_2/senanl_r0"))
# file.rename (from = paste0 (var$burden_estimate_folder, res_files),
#              to = paste0 ("previous_res/20230401/siareach_2/senanl_r0/", res_files))
