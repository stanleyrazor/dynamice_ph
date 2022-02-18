# plot_results.R
# update: 2022/02/17

library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(stringr)
library(readxl)
library(tidyr)

rm (list = c())
vac_stgs <- c("nomcv",             # (1) no vaccination
              "mcv1",              # (2) MCV1 only
              "mcv1-mcv2",         # (3) MCV1 + MCV2
              "mcv1-mcv2alt",      # (6) MCV1 + MCV2(alternative)
              "mcv1-sia",          # (5) MCV1 + SIA
              "mcv1-mcv2-sia",     # (4) MCV1 + MCV2 + SIA
              "mcv1-mcv2alt-sia",  # (7) MCV1 + MCV2(alternative) + SIA
              "mcv1-mcv2-siaplan"  # (8) MCV1 + MCV2 + SIA(plan)
)
vac_stg_names = c("No vaccination", "MCV1", "MCV1 + MCV2", "MCV1 + MCV2 (Alternative)",
                  "MCV1 + SIAs", "MCV1 + MCV2 + SIAs", "MCV1 + MCV2 (Alternative) + SIAs",
                  "MCV1 + MCV2 + SIAs (plan)")

eva_ctries = c("IND", "IDN", "NGA", "CHN", "PHL",
               "UGA", "ETH", "COD", "AGO", "NER",
               "PAK", "MDG", "SOM", "ZAF", "TZA",
               "MOZ", "TUR", "TCD", "BEN", "AFG")


# ------------------------------------------------------------------------------
## load and combine model outputs
# ------------------------------------------------------------------------------
# burden and vaccine doses
file_burden  <- NULL
for (scname in vac_stgs){
  scn_burden <- fread (paste0 (getwd(), "/central_burden_estimate/Portnoy/",
                               "central_burden_estimate_", scname, ".csv"))
  scn_burden [, `:=` (cases = cases0d + cases1d + cases2d,
                      deaths = deaths0d + deaths1d + deaths2d)]
  scn_burden [, comp := scname]
  file_burden <- rbind (file_burden, scn_burden)
}

cum_burden <- file_burden [, lapply (.SD, sum),
                           .SDcols = pops:deaths,
                           by = c("country", "year", "country_name", "comp")]


# check total cases by scenario
cum_burden [, .(total_cases = sum(cases)), by = "comp"]

# set up country name and ISO-3 code
country_names        <- cum_burden [year == 2000 & comp == "nomcv", country_name]
names(country_names) <- cum_burden [year == 2000 & comp == "nomcv", country]


# ------------------------------------------------------------------------------
## load WHO data for comparison
# ------------------------------------------------------------------------------
# WHO reported cases
input_WHOcase <- read_excel ("D:/research-data/incidence_series.xls", sheet = "Measles")
data_WHOcase  <- setDT (tidyr::pivot_longer(input_WHOcase, col = `2019`:`1980`,
                                            values_to = "cases", names_to = "year"))
data_WHOcase [,  `:=` (year = as.integer(year), comp = "nomcv")]
setnames (x = data_WHOcase,
          old = c("ISO_code", "Cname", "cases"),
          new = c("country", "country_name", "notifs"))

data_WHOcase [, country_name := country_names [country]]
data_WHOcase <- data_WHOcase [country_name %in% country_names & year >= 2000]


# ------------------------------------------------------------------------------
## plot overall cases and incidence trends
# ------------------------------------------------------------------------------
plt_measure <- function (sel_measure, sel_scns){
  plt <- ggplot (data = cum_burden [comp %in% sel_scns],
                 aes(x = year, y = get(sel_measure), colour = comp)) +
  geom_line (size = 0.9, alpha = 0.7) +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  labs (title = "Number of measles cases", x = "year", y = " ") +
  geom_point (data = data_WHOcase,
              aes(x = year, y = notifs, shape = "WHO reported cases"),
              size = 2, colour = "grey50") +
  scale_shape_manual ("", values = 4) +
  scale_colour_brewer("Scenarios", palette = "Dark2") +
  scale_y_log10 () +
  theme_bw () +
  theme (legend.position = "bottom") # legend.position = c(0.93, 0.75)
  return(plt)
}
ggsave ("plot/case-overview.pdf",
        plt_measure("cases", vac_stgs[c(1:3,5:6)]),
        height = 9, width = 15)

ggsave ("plot/case-overview-siaplan.pdf",
        plt_measure("cases", vac_stgs[c(6,8)]),
        height = 9, width = 15)

# ------------------------------------------------------------------------------
## plot % of zero-dose population at xth year of life
# ------------------------------------------------------------------------------
plt_zd_time <- function (sel_age, sel_ylife, sel_ctries, sel_palette){
  file_burden_y <- file_burden [age == sel_age & comp %in% vac_stgs[2:6]]
  file_burden_y [, `:=` (comp = factor (comp, levels = vac_stgs[2:6],
                                        labels = vac_stg_names[2:6]),
                         country_name = factor (country_name, levels = country_names[eva_ctries]))]
  plt <- ggplot (data = file_burden_y [country %in% sel_ctries],
                 aes(x = year, y = pops0d/pops, fill = comp)) +
    geom_col (position = "dodge") +
    scale_y_continuous (labels = scales::percent) +
    facet_wrap (vars(country_name), scales = "free", ncol = 1) +
    labs (title = paste0 ("% of zero-dose population in their ", sel_ylife, " year of life"),
          x = "year", y = " ") +
    scale_fill_brewer("Scenarios", palette = sel_palette) +
    theme_bw () +
    theme (legend.position = "bottom") # legend.position = c(0.93, 0.75)
  return(plt)
}


pdf (file = "plot/prc-zerodose-y3rd.pdf", width = 10, height = 10)
for (ig in 0:3){
  print (plt_zd_time (2, "third", eva_ctries [5*ig + (1:5)], "Blues"))
}
dev.off()

pdf (file = "plot/prc-zerodose-y2nd.pdf", width = 10, height = 10)
for (ig in 0:3){
  print (plt_zd_time (1, "second", eva_ctries [5*ig + (1:5)], "Greens"))
}
dev.off()


# ------------------------------------------------------------------------------
## plot dose distribution
# ------------------------------------------------------------------------------
pltdat_dose <- copy (cum_burden [, .(comp, country_name, country, year,
                                  doses, reachs0d, fvps)])
pltdat_dose [, `:=` (doses0 = reachs0d,
                     doses1 = fvps,
                     doses2 = doses-reachs0d-fvps)]
pltdat_dose <- setDT (pivot_longer (pltdat_dose [, !c("doses", "fvps", "reachs0d")],
                                    col = doses0:doses2,
                                    names_to = "measure", values_to = "value"))
pltdat_dose [, `:=` (country_name = factor (country_name,
                                            levels = country_names[eva_ctries]),
                     comp = factor (comp, levels = vac_stgs, labels = vac_stg_names),
                     measure = factor (measure,
                                       levels = c("doses0", "doses1", "doses2"),
                                       labels = c("0 dose", "1 dose", ">= 2 doses")))]

plt_dose <- function (sel_ctries){
  plt <- ggplot (data = pltdat_dose [country %in% sel_ctries & value > 0], # remove bar borders in the plot for zero values
                 aes(x = year, y = value/1e6)) +
    geom_col (aes(fill = measure), width = 1, colour = NA, size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_grid (rows = vars(country_name), cols = vars(comp), scales = "free_y") +
    labs (x = "year", y = "Number of doses (millions)") +
    scale_fill_brewer ("Population reached", palette = "Set2") +
    theme_bw ()
  return(plt)
}

pdf ("plot/dose.pdf", width = 14, height = 6)
for (ig in 0:3){
  print (plt_dose (eva_ctries[ig*5 + (1:5)]))
}
dev.off()

# sum over 2000-2020
pltdat_dose_sum <- pltdat_dose [, .(total_dose = sum(value)),
                                by = c("comp", "country_name", "country", "measure")]
plt_dose_sum <- function (sel_ctries, sel_scns){
  plt <- ggplot (data = pltdat_dose_sum [country %in% sel_ctries &
                                           comp %in% sel_scns &
                                           total_dose > 0], # remove bar borders in the plot for zero values
                 aes(x = comp, y = total_dose/1e6)) +
    geom_col (aes(fill = measure), width = 0.8, colour = NA, size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scale = "free_y", ncol = 5) +
    labs (x = "", y = "Number of doses over 2000-2020 (millions)") +
    scale_fill_brewer ("Population reached", palette = "Set2") +
    theme_bw () +

    theme (legend.position = "top",
           axis.text.x = element_text (angle = 60, hjust = 1))
  return(plt)
}
pdf ("plot/dose-sum.pdf", width = 13, height = 9)
print (plt_dose_sum (eva_ctries, vac_stg_names[c(2,3,5,6)]))
dev.off()


# ------------------------------------------------------------------------------
## plot cases and deaths by doses
# ------------------------------------------------------------------------------
# version 1: by time, scenario, country
plt_burden_byDose <- function (age_cutoff, sel_yrs, sel_scns, sel_ctries,
                               sel_title, sel_mea, change_to_pr = F){
  plt_data <- copy (file_burden [age >= age_cutoff & country %in% sel_ctries &
                                   year %in% sel_yrs])
  plt_data <- plt_data [, lapply (.SD, sum), .SDcols = cases0d:deaths2d,
                        by = c("country", "country_name", "comp", "year")]
  plt_data [, `:=` (prcases0d = cases0d/(cases0d+cases1d+cases2d),
                    prcases1d = cases1d/(cases0d+cases1d+cases2d),
                    prcases2d = cases2d/(cases0d+cases1d+cases2d),
                    prdeaths0d = deaths0d/(deaths0d+deaths1d+deaths2d),
                    prdeaths1d = deaths1d/(deaths0d+deaths1d+deaths2d),
                    prdeaths2d = deaths2d/(deaths0d+deaths1d+deaths2d))]
  plt_data <- setDT (pivot_longer (plt_data, col = cases0d:prdeaths2d,
                     names_to = "measure", values_to = "value"))
  plt_data [, `:=` (country_name = factor (country_name,
                                           levels = country_names[eva_ctries]),
                    comp = factor (comp, levels = vac_stgs, labels = vac_stg_names))]
  subset_data <- plt_data [str_sub(measure,1,5) == sel_mea & comp %in% sel_scns]
  subset_data [, `:=` (measure = factor (measure, labels = c("0 dose", "1 dose", ">=2 doses")))]

  plt <- ggplot (data = subset_data [value > 0], # remove bar borders in the plot for zero values
                 aes(x = year, y = value)) +
    geom_col (aes(fill = measure), width = 1, colour = NA,
              position = position_stack (reverse = TRUE)) +
    facet_grid (rows = vars(country_name), cols = vars(comp), scales = "free_y") +
    labs (x = "year", y = " ", title = sel_title) +
    scale_fill_brewer ("MCV received", palette = "Oranges") +
    theme_bw () +
    theme (legend.position = "right") # legend.position = c(0.93, 0.75)

  if (change_to_pr == T) {
    plt <- plt + scale_y_continuous (labels = scales::percent)
  }
  return(plt)
}

pdf ("plot/case-by-dose-age2up.pdf", width = 11, height = 6)
for (ic in 0:3){
  print (plt_burden_byDose (2, 2000:2020, vac_stg_names[c(2,3,5,6)],
                            eva_ctries[ic*5 + (1:5)],
                            "Number of cases, age >= 2", "cases", F))
}
dev.off()

# pdf ("plot/prcase-by-dose.pdf", width = 14, height = 6)
# for (ic in 0:3){
#   subset_data <- plt_data [country %in% eva_ctries[ic*5 + (1:5)]
#                            & str_sub(measure,1,7) == "prcases"]
#   print (plt_burden_byDose (subset_data, "Proportion of cases by doses", T))
# }
# dev.off()
#
# pdf ("plot/death-by-dose.pdf", width = 14, height = 6)
# for (ic in 0:3){
#   subset_data <- plt_data [country %in% eva_ctries[ic*5 + (1:5)]
#                            & str_sub(measure,1,6) == "deaths"]
#   print (plt_burden_byDose (subset_data, "Number of deaths", F))
# }
# dev.off()
#
# pdf ("plot/prdeath-by-dose.pdf", width = 14, height = 6)
# for (ic in 0:3){
#   subset_data <- plt_data [country %in% eva_ctries[ic*5 + (1:5)]
#                            & str_sub(measure,1,8) == "prdeaths"]
#   print (plt_burden_byDose (subset_data, "Proportion of deaths by doses", T))
# }
# dev.off()

# version 2: by scenario, country
plt_burden_byDose_sum <- function (age_cutoff, sel_yrs, sel_mea, sel_scns,
                                   sel_palette, sel_ylab){
  plt_data_sum <- copy (file_burden [age >= age_cutoff & year %in% sel_yrs])
  plt_data_sum <- plt_data_sum [, lapply (.SD, sum), .SDcols = cases0d:deaths2d,
                                by = c("country", "country_name", "comp")]
  plt_data_sum <- setDT (pivot_longer (plt_data_sum, col = cases0d:deaths2d,
                                       names_to = "measure", values_to = "value"))
  plt_data_sum [, `:=` (country_name = factor (country_name,
                                               levels = country_names[eva_ctries]),
                        comp = factor (comp, levels = vac_stgs, labels = vac_stg_names))]
  plt_data_sum [, measure_type := ifelse (str_sub(measure,1,4) == "case",
                                          "Cases", "Deaths")]
  plt_data_sum [, measure := factor (str_sub (measure,-2,-1),
                                     labels = c("0 dose", "1 dose", ">=2 doses"))]

  plt <- ggplot (data = plt_data_sum [measure_type == sel_mea &
                                        comp %in% sel_scns &
                                        value > 0], # remove bar borders in the plot for zero values
                 aes(x = comp, y = value/1e6, group = measure)) +
    geom_col (aes(fill = measure), width = 0.8, colour = NA, size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scales = "free_y", ncol = 5) +
    labs (x = "Scenario", y = sel_ylab) +
    scale_fill_brewer ("MCV received", palette = sel_palette) +
    theme_bw () +
    theme (legend.position = "top",
           axis.text.x = element_text (angle = 60, hjust = 1))
  return(plt)
}

pdf ("plot/case-by-dose-sum.pdf", width = 12, height = 12)
print (plt_burden_byDose_sum (0, 2000:2020, "Cases", vac_stg_names[c(2,3,5,6)],
                              "Oranges", "Cases over 2000-2020 (millions)"))
dev.off()

pdf ("plot/case-by-dose-sum-2016to2020-age2up.pdf", width = 13, height = 9)
print (plt_burden_byDose_sum (2, 2016:2020, "Cases", vac_stg_names[c(2,3,5,6)],
                              "Oranges", "Cases over 2016-2020 (millions), age 2+"))
dev.off()

# pdf ("plot/death-by-dose-sum.pdf", width = 12, height = 12)
# print (plt_burden_byDose_sum (0, "Deaths", vac_stg_names[c(2,3,5,6)], "Purples"))
# dev.off()
