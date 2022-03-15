# process_coverage.R
# generate vaccine coverage files
# update: 2022/03/01

library(data.table)
library(readxl)
library(stringr)
library(ggplot2)
library(scales)

rm (list = c())

# top 20 countries with measles cases over 2017-2019
# covering 83% of the total global cases
# India, Indonesia, Nigeria, China, Philippines,
# Uganda, Ethiopia, Democratic Republic of the Congo, Angola, Niger,
# Pakistan, Madagascar, Somalia, South Africa, United Republic of Tanzania,
# Mozambique, Turkey (not included in VIMC), Chad, Benin, Afghanistan
sel_ctries = c("IND", "IDN", "NGA", "CHN", "PHL",
               "ETH", "UGA", "AGO", "COD", "MOZ",
               "SOM", "PAK", "ZAF", "MDG", "NER",
               "TZA", "TUR", "TCD", "BEN", "AFG")

# set up column names for the output file
# based on the format of VIMC coverage file
# "D:/dynamice_private/vaccine_coverage/coverage_202110gavi-3_measles-campaign-default.csv"
sel_cols <- c("vaccine", "dose", "activity_type", "country_code", "country", "year",
              "age_first", "age_last", "age_range_verbatim", "target", "coverage")


# ------------------------------------------------------------------------------
## routine immunisation (MCV1, MCV2)
# ------------------------------------------------------------------------------
# load WHO-UNICEF data MCV1 and MCV2
# https://immunizationdata.who.int/pages/coverage/mcv.html

cov_who_mcv <- read_excel ("D:/research-data/Measles vaccination coverage.xlsx", sheet = 1)
cov_who_mcv <- setDT(cov_who_mcv)[(CODE %in% sel_ctries) & (COVERAGE_CATEGORY == "WUENIC")]
cov_who_mcv [is.na(COVERAGE), COVERAGE := 0]

# merge with template file
cov_file_routine <- cov_who_mcv [, .SD, .SDcols = c("CODE", "NAME", "YEAR", "ANTIGEN", "COVERAGE")]
setnames (x = cov_file_routine, new = c("country_code", "country", "year", "vaccine", "coverage"))
cov_file_routine [, `:=` (coverage           = coverage/100,
                          activity_type      = "routine",
                          age_first          = numeric(),
                          age_last           = numeric(),
                          age_range_verbatim = character(),
                          target             = numeric())]
cov_file_routine [, dose := vaccine]


# add data for missing years
sim_years <- 1980:2020
for (ivac in c("MCV1","MCV2")){
  for(ictry in sel_ctries){
    single_vac_ctry  <- cov_file_routine [vaccine == ivac & country_code == ictry]
    miss_yrs <- sim_years[!(sim_years %in% single_vac_ctry$year)]

    if (length (miss_yrs) > 0){
      for (imiss in miss_yrs){
         cov_file_routine <- rbind (cov_file_routine, copy(single_vac_ctry[1][, `:=` (year = imiss, coverage = 0)]))
      }
    }
    remove(single_vac_ctry, miss_yrs)
  }
}
setorder (cov_file_routine, vaccine, country_code, year)

# # get expected coverage without COVID-related disruptions in 2020
# # by assigning 2020 estimate or 2019 estimate +1%, whichever is larger
# cov_file_routine_2020 <- cov_file_routine [year == 2020]
# cov_file_routine_2020 <- cov_file_routine_2020 [cov_file_routine [year == 2019],
#                                                 .(country_code, year, vaccine, coverage, i.coverage),
#                                                 on = .(vaccine = vaccine , country_code = country_code)]
# cov_file_routine_2020 [, coverage_adj := ifelse ( (coverage-i.coverage) > 0.01, coverage, (i.coverage+0.01))]
# cov_file_routine [cov_file_routine_2020 , coverage := coverage_adj,
#                   on = .(vaccine = vaccine, country_code = country_code, year = year)]
# cov_file_routine [coverage > 0.95, coverage := 0.95]
#
# # get future projections
# routine_type <- c("MCV1", "MCV2")
# for (vaccine_type in routine_type){
#   cov_future_routine <- rbindlist (lapply (sel_ctry, function (ictry) {
#     rbindlist (lapply (2021:2100, function (iyr) {
#       copy (cov_file_routine [vaccine == vaccine_type & country_code == ictry & year == 2020, ])[
#         , `:=`(coverage = min (0.95, coverage + 0.01*(iyr-2020)), year = iyr)]}))
#   }))
#   cov_file_routine <- rbind (cov_file_routine, cov_future_routine)
# }


# ------------------------------------------------------------------------------
## supplementary immunisation activities (SIAs)
# ------------------------------------------------------------------------------
# load WHO data SIAs
# https://www.who.int/entity/immunization/monitoring_surveillance/data/Summary_Measles_SIAs.xls
cov_who_sia <- as.data.table (read_excel ("D:/research-data/Summary_Measles_SIAs.xls",
                                          sheet = 2, skip = 1)) # "SIAs_Jan2000_to_Dec2020"

# exclude campaigns that were not carried out
table (cov_who_sia [ISO %in% sel_ctries, `Implementation status`])
# done   planned postponed
# 375         5         9
cov_who_sia <- cov_who_sia [ISO %in% sel_ctries & `Implementation status` == "done"]

# ensure consistent format
cov_who_sia [, `:=` (`Age group` = tolower (gsub("[[:space:]]",
                                                 "", `Age group`)),
                     Extent = tolower (gsub("[[:space:]]", "", Extent)))]

# ------------------------------------------------------------------------------
# # summarise basic information
# table (cov_who_sia$Activity)
# # Campaign    CatchUp     Child health days       FollowUp        HRA
# # 20          117         4                       159             1
# # MopUp       NIDs        Outbreak Response       SIA
# # 26          1           46                      1
# table (cov_who_sia$Extent)
# # national rollover-nat sub-national      unknown
# # 85           87          199            4
# table (cov_who_sia$Intervention)
# # Measles   MMR     MR
# # 316       13      46
# table (cov_who_sia$Country)
# # Afghanistan  Angola     Benin         Chad         China         DRCongo
# # 18           9          8             16           13            59
# # Ethiopia     India      Indonesia     Madagascar   Mozambique    Niger
# # 25           41         24            9            7             7
# # Nigeria      Pakistan   Philippines   Somalia      South Africa  Tanzania
# # 16           31         18            31           11            10
# # Turkey       Uganda
# # 12           10
# ------------------------------------------------------------------------------

# rename columns
setnames (cov_who_sia,
          old = c("Year", "Country", "ISO", "Age group", "Target population",
                  "Reached population", '% Reached'),
          new = c("year", "country", "country_code", "age_range_verbatim", "target",
                  "reach_pop", "coverage"))

# adjust format for following extraction of precise age
table (cov_who_sia$age_range_verbatim)
cov_who_sia [age_range_verbatim == "<5y",
             age_range_verbatim := "6m-4y"]
cov_who_sia [age_range_verbatim == ">12yoroutsidetargetvaccinated",
             age_range_verbatim := "13y-100y"]
cov_who_sia [age_range_verbatim == "school-age",
             age_range_verbatim := "6y-17y"]
cov_who_sia [age_range_verbatim %in% c("airportworkers", "tourismworkers"),
             age_range_verbatim := "15y-65y"] # general age range for labour force
cov_who_sia [age_range_verbatim == "refugees",
             age_range_verbatim := "6m-100y"] # all eligible age groups
cov_who_sia [age_range_verbatim %in% c("military(1980-91cohorts)", "hcw(1980-91cohorts)"),
             age_range_verbatim := "29y-40y"] # calendar year 2019
cov_who_sia [age_range_verbatim == "unknown",
             age_range_verbatim := "9-59m"]   # Chad, inferred from historical data
cov_who_sia [country_code == "ZAF" & `Areas/comments` == "5y - <20yrs in Tshwane District",
             age_range_verbatim := "5y-19y"]
cov_who_sia [country_code == "ZAF" & `Areas/comments` == "9m - <20yrs in 5 districts in Gauteng Province",
             age_range_verbatim := "9m-19y"]  # based on comments
cov_who_sia [country_code == "COD" & `Areas/comments` == "or 6 M-9 Y",
             age_range_verbatim := "6m-9y"]

cov_who_sia [, `:=`(age_first_verb = stringr::str_extract (age_range_verbatim, "[^-]+"),
                    age_last_verb  = stringr::str_extract (age_range_verbatim, "[^-]+$"))]

cov_who_sia [, `:=`(age_first_num  = as.double (stringr::str_extract (age_first_verb, "\\d+")),
                    age_first_unit = stringr::str_extract (age_first_verb, "\\D+$"),
                    age_last_num   = as.double (stringr::str_extract (age_last_verb, "\\d+")),
                    age_last_unit  = stringr::str_extract (age_last_verb, "\\D+$"))]

cov_who_sia [stringr::str_detect (age_range_verbatim, "<"), age_last_num := age_last_num - 1]
cov_who_sia [is.na(age_first_unit), age_first_unit := age_last_unit]

# use precise age for those <=3 years old for weekly age structure
# treat 36-month old as the end of 2-year old
# age definition used in UNWPP data: X <= age < (X+1), X = 0, 1, ...100
trans_m_to_y <- function (mage){
  yage <- ifelse (mage < 36, mage/12, ifelse (mage == 36, 2.999999, floor(mage/12)))
  return (yage)}
cov_who_sia [age_first_unit == "m", age_first := trans_m_to_y(age_first_num)]
cov_who_sia [age_last_unit  == "m", age_last  := trans_m_to_y(age_last_num)]
cov_who_sia [age_first_unit == "y", age_first := age_first_num]
cov_who_sia [age_last_unit  == "y", age_last  := age_last_num]

# adjust age ranges with operators (usually occur in upper bounds)
cov_who_sia [stringr::str_detect (age_range_verbatim, "<") & age_last_unit == "y", age_last := age_last -1]
cov_who_sia [stringr::str_detect (age_range_verbatim, ">") & age_last_unit == "y", age_last := 100]

# deliver SIA to children older than 6 months old, as WHO recommended
cov_who_sia [age_first < 0.5, age_first := 0.5]
cov_who_sia [age_last < 0.5, age_last := 0.5]

# check age ranges for SIAs
cov_who_sia [cov_who_sia[,.I[1], by = age_range_verbatim]$V1, age_first_verb:age_last]


## calculate the size of target population
# an age group of x in the UNWPP data covers population aged [x, x+1)
load (file = "data/data_pop.rda")
get_targetpop <- function (icty, iyr, iage1, iage2){
  if ((iage1 < 1) & (iage2 < 1)) {
    targetpop <- (iage2-iage1)*data_pop [country_code == icty &
                                           year == as.integer(iyr) &
                                           age_from == 0, value]
  } else {
    targetpop <- sum (data_pop [country_code == icty & year == as.integer(iyr) &
                                  age_from %in% (ceiling(iage1):floor(iage2)), value])

    if (ceiling(iage1) > iage1){
      targetpop <- targetpop +
        (ceiling(iage1)-iage1)*data_pop [country_code == icty & year == as.integer(iyr) &
                                           age_from == floor(iage1), value] }

    if (iage2 > floor(iage2)){
      targetpop <- targetpop -
        (ceiling(iage2)-iage2)*data_pop [country_code == icty & year == as.integer(iyr) &
                                           age_from == floor(iage2), value]}
  }
  return (targetpop)
}


# ## replace missing data based on historical data - check individually
# # infer the reached population in the same country and for similar type of
# # activity and age group if the data is available
# # basic rules: target population among total population, reached population
# # among total population, or the administration coverage remain unchanged
# cov_who_sia [is.na(reach_pop)]
#
# # 2019, DR Congo, FollowUp -> inferred from 2017, DR Congo, FollowUp
# cov_who_sia [country_code  == "COD" & Activity == "FollowUp" & year == 2019,
#              reach_pop := 1.035 * target]
#
# # 2020, DR Congo, Outbreak response -> inferred from 2015, DR Congo, Outbreak response, 6m-10y
# cov_who_sia [country_code  == "COD" & Activity == "Outbreak Response" & year == 2020,
#              reach_pop := 1.00 * target]
#
# # 2020, Tanzania, Outbreak response -> inferred from 2020, DR Congo, Outbreak response
# # 1268185/get_targetpop("COD", 2020, 0.5, 9) = 0.0455901
# cov_who_sia [country_code == "TZA" & Activity == "Outbreak Response" & year == 2020,
#              target := 0.0455901*get_targetpop(country_code, year, age_first, age_last)]
# cov_who_sia [country_code == "TZA" & Activity == "Outbreak Response" & year == 2020,
#              reach_pop := 1.00 * target]
#
# # 2002, Somalia, MopUp -> inferred from 2011, Somalia, MopUp, 9-59m
# # 2924/get_targetpop("SOM", 2011, 0.75, 4) = 0.001498303
# cov_who_sia [country_code == "SOM" & Activity == "MopUp" & year == 2002,
#              reach_pop := 0.001498303*get_targetpop(country_code, year, age_first, age_last)]
#
# # 2020, Somalia, FollowUp -> inferred from 2019, Somalia, FollowUp
# cov_who_sia [country_code == "SOM" & Activity == "FollowUp" & year == 2020,
#              reach_pop := 0.92*target]
#
# # 2003, Turkey, CatchUp, 9m-14y -> inferred from 2003, Turkey, CatchUp, 6-14y
# cov_who_sia [country_code == "TUR" & Activity == "CatchUp" & year == 2003,
#              reach_pop := 0.97*target]
#
# # 2001 & 2005, Philippines, FollowUp -> inferred from 2002, Philippines, FollowUp
# # 507463/get_targetpop("PHL", 2002, 0.75, 4) = 0.05415916
# cov_who_sia [country_code == "PHL" & Activity == "FollowUp" & year == 2001,
#              target := 0.05415916*get_targetpop(country_code, year, age_first, age_last)]
# cov_who_sia [country_code == "PHL" & Activity == "FollowUp" & year == 2001,
#              reach_pop := 0.99*target]
# cov_who_sia [country_code == "PHL" & Activity == "FollowUp" & year == 2005,
#              target := 0.05415916*get_targetpop(country_code, year, age_first, age_last)]
# cov_who_sia [country_code == "PHL" & Activity == "FollowUp" & year == 2005,
#              reach_pop := 0.99*target]
#
# # 2013, Philippines, Outbreak Response -> inferred from 2010, Philippines, High risk areas (HRA)
# # 420129/get_targetpop("PHL", 2010, 0.75, 4) = 0.04507288
# cov_who_sia [country_code == "PHL" & Activity == "Outbreak Response" & year == 2013,
#              reach_pop := 0.04507288*get_targetpop(country_code, year, age_first, age_last)]

# remove missing data
cov_who_sia <- cov_who_sia [!is.na(reach_pop)]
dim(cov_who_sia)
# 366  29

## calculate target population
cov_who_sia [, target := get_targetpop (country_code, year, age_first, age_last),
             by = seq_len (nrow(cov_who_sia))]

# calculate coverage
cov_who_sia [, coverage := reach_pop/target]
cov_who_sia [coverage > 1, coverage := 1]

summary(cov_who_sia$coverage)
# Min.   1st Qu.    Median      Mean      3rd Qu.   Max.
# 0.0000319 0.0257182 0.1374162 0.3145602 0.5425627 1.0000000


## plot general trend
plt_sia <- ggplot (data = cov_who_sia,
                   aes (x = year, y = coverage, colour = Activity, shape = Extent)) +
  scale_x_continuous (breaks = pretty_breaks ()) +
  facet_wrap(vars(country), ncol = 4) +
  geom_point (alpha = 0.75, size = 1.2) + #, fill = NA
  scale_shape_manual (values = c(0:3)) +
  labs (title = "SIA coverage", x = "Year", y = "Vaccine coverage") +
  theme_bw () +
  theme (axis.text.x = element_text(size = 8),
         axis.text.y = element_text(size = 8))
ggsave ("plot/coverage-sia.pdf", plt_sia, width = 28, height = 16, units = "cm")

# # get future projections - every three years
# for (ictry in sel_ctry){
#   year_last_sia    <- max (cov_sia_full [ISO == ictry, year])
#   years_future_sia <- seq (year_last_sia, 2100, 3)
#
#   cov_future_sia <- rbindlist (lapply (years_future_sia[-1], function (iyr){
#     newyr_row <- copy (cov_sia_full [ISO == ictry & year == year_last_sia])[ , year := iyr]
#     newyr_row [, target := get_targetpop (ISO, year, age_first, age_last),
#                by = seq_len (nrow(newyr_row))]
#     return (newyr_row)
#   }))
#   cov_sia_full <- rbind (cov_sia_full, cov_future_sia)
# }

# extract the month of SIA implementation
cov_who_sia [, `:=` (start_m = as.numeric (format(`Start date`,"%m")))]
# cov_who_sia [is.na (start_m)][, .N, by = c("country", "year")]
# assume missing values = mid-year
cov_who_sia [is.na (start_m), start_m := 6.5]

# merge with template coverage file
setnames (x = cov_who_sia, old = c("Intervention"), new = c("vaccine"))
cov_who_sia [, `:=` (dose = vaccine, activity_type = "campaign")]
cov_who_sia [, vaccine := "SIA"]
cov_file_sia <- cov_who_sia [, .SD, .SDcols = c(sel_cols, "start_m")]
setorder (cov_file_sia, vaccine, country_code, year)

# add a version with planned SIAs only (no SIAs for outbreak response)
cov_who_sia_plan  <- cov_who_sia [Activity != "Outbreak Response"]
cov_file_sia_plan <- cov_who_sia_plan [, .SD, .SDcols = c(sel_cols, "start_m")]
setorder (cov_file_sia_plan, vaccine, country_code, year)
dim(cov_file_sia_plan)


# ------------------------------------------------------------------------------
## output coverage files
# ------------------------------------------------------------------------------
outfile_mcv1_mcv2_sia <- rbind (cov_file_routine, cov_file_sia, fill = TRUE)
outfile_mcv1_mcv2_sia [country == "DRCongo", country := "Democratic Republic of the Congo"]
outfile_mcv1_mcv2_sia [country == "Tanzania", country := "United Republic of Tanzania"]
outfile_mcv1_mcv2_sia [, scenario := "mcv1-mcv2-sia"]
fwrite (outfile_mcv1_mcv2_sia, "coverage/coverage_mcv1-mcv2-sia.csv")

outfile_mcv1_mcv2 <- copy (outfile_mcv1_mcv2_sia) [activity_type == "routine"]
outfile_mcv1_mcv2 [, scenario := "mcv1-mcv2"]
fwrite (outfile_mcv1_mcv2, "coverage/coverage_mcv1-mcv2.csv")

outfile_mcv1 <- copy (outfile_mcv1_mcv2) [vaccine == "MCV2", coverage := 0]
outfile_mcv1 [, scenario := "mcv1"]
fwrite (outfile_mcv1, "coverage/coverage_mcv1.csv")

outfile_mcv1_sia <- copy (outfile_mcv1_mcv2_sia) [vaccine == "MCV2", coverage := 0]
outfile_mcv1_sia [, scenario := "mcv1-sia"]
fwrite (outfile_mcv1_sia, "coverage/coverage_mcv1-sia.csv")

outfile_nomcv <- copy (outfile_mcv1_mcv2) [, coverage := 0]
outfile_nomcv [, scenario := "nomcv"]
fwrite (outfile_nomcv, "coverage/coverage_nomcv.csv")

# check years and months of implementation
temp_sum <- outfile_mcv1_mcv2_sia [activity_type == "campaign"][, .N, by = c("year", "start_m", "country_code")]
temp_sum [N > 1]

# add alternative SIA data
outfile_mcv1_mcv2_sia_plan <- rbind (cov_file_routine, cov_file_sia_plan, fill = TRUE)
outfile_mcv1_mcv2_sia_plan [country == "DRCongo", country := "Democratic Republic of the Congo"]
outfile_mcv1_mcv2_sia_plan [country == "Tanzania", country := "United Republic of Tanzania"]
outfile_mcv1_mcv2_sia_plan [, scenario := "mcv1-mcv2-siaplan"]
fwrite (outfile_mcv1_mcv2_sia_plan, "coverage/coverage_mcv1-mcv2-siaplan.csv")


# ------------------------------------------------------------------------------
## add alternative MCV2 scenario - 1
# ------------------------------------------------------------------------------
outfile_mcv1_mcv2 <- fread ("coverage/coverage_mcv1-mcv2.csv")
dat_mcv2_2020 <- outfile_mcv1_mcv2 [vaccine == "MCV2" & year == 2020]
dat_mcv1_2020 <- outfile_mcv1_mcv2 [vaccine == "MCV1" & year == 2020]
dat_mcv2alt  <- dat_mcv2_2020 [dat_mcv1_2020 [, .(country_code, coverage)],
                               on = .(country_code = country_code)]
dat_mcv2alt [coverage == 0, coverage := i.coverage - 0.1]

# 'mcv1-mcv2alt' scenario
outfile_mcv1_mcv2alt <- copy(outfile_mcv1_mcv2) [dat_mcv2alt [, .(country_code, coverage)],
                                                 on = .(country_code = country_code)]
outfile_mcv1_mcv2alt [vaccine == "MCV2", coverage := ifelse (year < 2000, 0, i.coverage)]
outfile_mcv1_mcv2alt [, i.coverage := NULL]
fwrite (outfile_mcv1_mcv2alt, "coverage/coverage_mcv1-mcv2alt.csv")

# 'mcv1-mcv2alt-sia' scenario
outfile_mcv1_mcv2_sia <- fread ("coverage/coverage_mcv1-mcv2-sia_old.csv")
outfile_mcv1_mcv2alt_sia  <- copy(outfile_mcv1_mcv2_sia) [dat_mcv2alt [, .(country_code, coverage)],
                                                 on = .(country_code = country_code)]
outfile_mcv1_mcv2alt_sia [vaccine == "MCV2", coverage := ifelse (year < 2000, 0, i.coverage)]
outfile_mcv1_mcv2alt_sia [, i.coverage := NULL]
fwrite (outfile_mcv1_mcv2alt_sia, "coverage/coverage_mcv1-mcv2alt-sia_old.csv")


# ------------------------------------------------------------------------------
## add alternative MCV2 scenario - 2
# ------------------------------------------------------------------------------
outfile_mcv1_mcv2 <- fread ("coverage/coverage_mcv1-mcv2.csv")
dat_mcv2_yrs <- outfile_mcv1_mcv2 [vaccine == "MCV2"]
dat_mcv1_yrs <- outfile_mcv1_mcv2 [vaccine == "MCV1"]
dat_mcv2alt_yrs  <- dat_mcv2_yrs [dat_mcv1_yrs [, .(country_code, year, coverage)],
                                  on = .(country_code = country_code,
                                         year = year)]
dat_mcv2alt_yrs [year >= 2000, coverage := i.coverage - 0.1]

# 'mcv1-mcv2alt' scenario
outfile_mcv1_mcv2alt <- copy(outfile_mcv1_mcv2) [dat_mcv2alt_yrs [, .(country_code, year, coverage)],
                                                 on = .(country_code = country_code,
                                                        year = year)]
outfile_mcv1_mcv2alt [vaccine == "MCV2", coverage := ifelse (coverage>=i.coverage, coverage, i.coverage)]
outfile_mcv1_mcv2alt [, i.coverage := NULL]
fwrite (outfile_mcv1_mcv2alt, "coverage/coverage_mcv1-mcv2alt.csv")

# 'mcv1-mcv2alt-sia' scenario
outfile_mcv1_mcv2_sia <- fread ("coverage/coverage_mcv1-mcv2-sia.csv")

outfile_mcv1_mcv2alt_sia  <- rbind (copy (outfile_mcv1_mcv2_sia [vaccine != "MCV2"]),
                                    copy (outfile_mcv1_mcv2alt [vaccine == "MCV2"]))
fwrite (outfile_mcv1_mcv2alt_sia, "coverage/coverage_mcv1-mcv2alt-sia.csv")


# ------------------------------------------------------------------------------
## plot general trend
# ------------------------------------------------------------------------------
outfile_mcv1_mcv2_sia <- fread ("coverage/coverage_mcv1-mcv2-sia.csv")
outfile_mcv1_mcv2alt  <- fread ("coverage/coverage_mcv1-mcv2alt.csv")
plt_data <- rbind (outfile_mcv1_mcv2_sia,
                   copy (outfile_mcv1_mcv2alt [vaccine == "MCV2"])[, vaccine := "MCV2 (alternative)"])

# update country names
plt_data [country_code == "COD", country := "DRC"]
plt_data [country_code == "TZA", country := "Tanzania"]

# rank countries by IHME burden
country_names        <- plt_data [year == 2000 & vaccine == "MCV1", country]
names(country_names) <- plt_data [year == 2000 & vaccine == "MCV1", country_code]
plt_data [, country := factor (country, levels = country_names[sel_ctries])]

pdf ("plot/coverage-check.pdf", width = 12, height = 8)
plt_cov <- ggplot (data = plt_data [vaccine != "SIA"],
                   aes (x = year, y = coverage, colour = vaccine, linetype = vaccine)) +
  scale_x_continuous (breaks = pretty_breaks ()) +
  geom_line (size = 1) +
  facet_wrap (vars(country), nrow = 4) +
  labs (title = " ", x = "Year", y = "Vaccine coverage") +
  theme_bw() +
  theme (legend.position  = "bottom",
         legend.direction = "horizontal",
         legend.key.size = unit (1.2, 'cm'),
         legend.text = element_text (size = 10),
         axis.text.x = element_text (angle = 60, hjust = 1),
         strip.text.x = element_text (size = 10),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         plot.margin = margin (0.1, 0.1, 0, 0.2, "cm"))
plt_cov <- plt_cov +
  geom_point (data = plt_data [vaccine == "SIA"], aes (x = year, y = coverage)) +
  scale_colour_manual ("Dose", #c("#253582ff", "#b8627dff", "#b8627dff", "#f9b641ff")
                       values = c("#42b540", "#00468b", "#0099b4", "#ed0000")) +
  scale_linetype_manual ("Dose", values = c(1,1,2,0)) +
  guides(color = guide_legend (override.aes = list (shape = c(NA,NA,NA,16))))
print(plt_cov)
dev.off()

