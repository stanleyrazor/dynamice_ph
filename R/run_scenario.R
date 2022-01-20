# run_scenario
# execute simulations for MR-MAPs scenarios
# update: 2022/01/20

rm(list = ls())

# load libraries, data, functions
library (data.table)
library (stringr)
library (ggplot2)
library (scales)
library (foreach)
library (iterators)
library (parallel)
library (doParallel) #Loading required package: iterators
library (countrycode)
library (ggpubr)

# rcpp codes
library (Rcpp)
sourceCpp("Rcpp/rcpp_spinup.cpp")
sourceCpp("Rcpp/rcpp_vaccine_oney.cpp")

load(file = "data/data_pop.rda")
load(file = "data/data_cfr_portnoy_21.rda")
load(file = "data/data_contact_syn.rda")
load(file = "data/data_r0.rda")      # assume fixed R0
load(file = "data/data_timeliness_maps.rda")  # to include Turkey
load(file = "data/data_lexp_remain_maps.rda") # to include Turkey
load(file = "data/data_template.rda")

# assume a fixed R0 for the central run
# median R0 for least developed countries, vaccine era, from Guerra et al. (2017)
adj.fixR0 = 15.9  # NA
if (!is.na(adj.fixR0)){
  data_r0 [ , r0 := adj.fixR0]
  # add r0 for Turkey
  data_r0 <- rbind (data_r0,
                    copy(data_r0[1])[, `:=` (country = "Turkey",
                                             country_code = "TUR")])
}

# assume Turkey has the same CFR as Syria
# use Portnoy 2021 estimates - no difference by scenario
data_cfr_portnoy_21 <- rbind (data_cfr_portnoy_21,
                              copy(data_cfr_portnoy_21[country == "SYR"])[, country := "TUR"])

source ("R/logs.R")
source ("R/utils.R")
source ("R/functions_rcpp.R")

# set up variables
var <- list (
  # vaccine coverage
  vaccine_coverage_folder           = "coverage/",
  coverage_prefix                   = "coverage",
  touchstone                        = "_",
  antigen                           = NULL,
  vaccine_coverage_subfolder        = "scenarios/",

  # disease burden
  burden_estimate_folder            = "central_burden_estimate/",

  # diagnostic plots folder
  plot_folder                       = "plot/",

  # modelling group name
  group_name                        = NULL,

  # log file name
  log_name                          = "test_log",

  # countries - specify iso3 codes to analyse only these countries
  #             or set it to "all" to analyse all included countries
  countries                         =  c("IND", "IDN", "NGA", "CHN", "PHL",
                                         "UGA", "ETH", "COD", "AGO", "NER",
                                         "PAK", "MDG", "SOM", "ZAF", "TZA",
                                         "MOZ", "TUR", "TCD", "BEN", "AFG"),
  cluster_cores                     = 1,    # number of cores
  psa                               = 0     # psa runs; 0 for single central run
)

dir.create (file.path (paste0 (getwd(), "/", var$burden_estimate_folder, "Portnoy/")), recursive = T)


# prepare coverage inputs
vaccine_strategies <- c("nomcv",         # (1) no vaccination
                        "mcv1",          # (2) MCV1 only
                        "mcv1-mcv2",     # (3) MCV1 + MCV2
                        "mcv1-mcv2-sia", # (4) MCV1 + MCV2 + SIA
                        "mcv1-sia",      # (5) MCV1 + SIA
                        "mcv1-mcv2alt"   # (6) MCV1 + MCV2(alternative)
)

# set SIAs and vaccination parameters for each scenario to minimize errors for running
set_sia         <- c (0, 0, 0, 1, 1, 0)
set_vaccination <- c (0, 1, 2, 2, 1, 2)

# prepare coverage input data - update when the data are changed
adj.covfiles <- FALSE
if(adj.covfiles){
  # Generate 2 coverage input files for routine and SIA vaccination
  for (index in 1:length(vaccine_strategies)) {
    create_vaccine_coverage_routine_sia (
      vaccine_coverage_folder    = var$vaccine_coverage_folder,
      vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
      coverage_prefix            = var$coverage_prefix,
      touchstone                 = var$touchstone,
      antigen                    = var$antigen,
      scenario_name              = vaccine_strategies [index],
      rev_cov                    = FALSE
    )
  }
}

# run model
for (index in 1:length(vaccine_strategies)){

  scenario_name  <- vaccine_strategies [index]
  print (scenario_name)
  scenario_number <- sprintf ("scenario%02d", index)

  # run model and estimate cases
  burden_estimate_file <- runScenario_rcpp (
    vaccine_coverage_folder    = var$vaccine_coverage_folder,
    vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
    coverage_prefix            = var$coverage_prefix,
    antigen                    = var$antigen,
    scenario_name              = scenario_name,
    save_scenario              = scenario_number,
    burden_estimate_folder     = var$burden_estimate_folder,
    group_name                 = var$group_name,
    log_name                   = var$log_name,
    countries                  = var$countries,
    cluster_cores              = var$cluster_cores,
    psa                        = var$psa,
    vaccination                = set_vaccination [index],
    using_sia                  = set_sia         [index],
    contact_mat                = "syn",
    sim_years                  = 1980:2020
  )

  # separately estimate dalys
  burden_estimate_file <- paste0 ("central_burden_estimate_",
                                  var$antigen,
                                  var$group_name,
                                  scenario_name,
                                  ".csv")

  # merge outputs into csv files
  merge_case_csv (vaccine_coverage_folder    = var$vaccine_coverage_folder,
                  vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
                  antigen                    = var$antigen,
                  scenario_name              = vaccine_strategies [index],
                  save_scenario              = sprintf ("scenario%02d", index),
                  burden_estimate_folder     = var$burden_estimate_folder,
                  group_name                 = var$group_name,
                  log_name                   = var$log_name,
                  psa                        = var$psa,
                  vaccination                = set_vaccination [index],
                  using_sia                  = set_sia         [index],
                  folder_date                = "20220120",
                  sim_years                  = 1980:2020,
                  chunksize                  = 1)

  # estimate deaths and DALYs by Portnoy's methods
  estimate_deaths_dalys_21 (antigen                = var$antigen,
                            group_name             = var$group_name,
                            scenario_name          = vaccine_strategies  [index],
                            log_name               = var$log_name,
                            burden_estimate_folder = var$burden_estimate_folder,
                            psa                    = var$psa,
                            chunksize              = 1)
}
