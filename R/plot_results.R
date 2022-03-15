# plot_results.R
# update: 2022/03/09

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
               "ETH", "UGA", "AGO", "COD", "MOZ",
               "SOM", "PAK", "ZAF", "MDG", "NER",
               "TZA", "TUR", "TCD", "BEN", "AFG")


# ------------------------------------------------------------------------------
## load and combine model outputs
# ------------------------------------------------------------------------------
# burden and vaccine doses
file_burden  <- NULL
for (scname in vac_stgs[1:7]){
  scn_burden <- fread (paste0 (getwd(), "/central_burden_estimate/Portnoy/",
                               "central_burden_estimate_", scname, ".csv"))
  scn_burden [, `:=` (cases = cases0d + cases1d + cases2d,
                      deaths = deaths0d + deaths1d + deaths2d)]
  scn_burden [, comp := scname]
  file_burden <- rbind (file_burden, scn_burden)
}

file_burden [country == "COD", country_name := "DRC"]
file_burden [country == "TZA", country_name := "Tanzania"]

cum_burden <- file_burden [, lapply (.SD, sum),
                           .SDcols = pops:deaths,
                           by = c("country", "year", "country_name", "comp")]

# set up country name and ISO-3 code
country_names        <- cum_burden [year == 2000 & comp == "nomcv", country_name]
names(country_names) <- cum_burden [year == 2000 & comp == "nomcv", country]

# calculate incidence rate (case per million)
cum_burden [, incrateM := (cases/pops)*1e6]


# ------------------------------------------------------------------------------
## load WHO data for comparison
# ------------------------------------------------------------------------------
# WHO reported cases
input_WHOcase <- read_excel ("D:/research-data/incidence_series.xls", sheet = "Measles")
data_WHOcase  <- setDT (tidyr::pivot_longer(input_WHOcase, col = `2019`:`1980`,
                                            values_to = "cases", names_to = "year"))
data_WHOcase [, year := as.integer(year)]
setnames (x = data_WHOcase,
          old = c("ISO_code", "Cname", "cases"),
          new = c("country", "country_name", "notifs"))
data_WHOcase <- data_WHOcase [cum_burden [comp == "nomcv" & year <= 2019,
                                          .(comp, country, year, pops)],
                              on = .(country = country, year = year)]
data_WHOcase [, incrateM := 1e6*(notifs/pops)]

data_WHOcase [, country_name := country_names [country]]
data_WHOcase <- data_WHOcase [country_name %in% country_names & year >= 2000]


# ------------------------------------------------------------------------------
## calculate averted burden and number needed to treat (NNV)
# ------------------------------------------------------------------------------
total_burden <- cum_burden [, .(total_cases = sum(cases),
                                total_deaths = sum(deaths),
                                total_dalys = sum(dalys),
                                total_doses = sum(doses)),
                            by = c("comp", "country_name", "country")]
# add summary estimates of 20 countries
total_burden <- rbind (total_burden,
                       total_burden [,  .(country_name = "Global",
                                          country = "Global",
                                          total_cases = sum(total_cases),
                                          total_deaths = sum(total_deaths),
                                          total_dalys = sum(total_dalys),
                                          total_doses = sum(total_doses)),
                                          by = c("comp")])
total_burden [country == "Global"]

# get averted burden, additional doses, and NNV
cal_avtnnv <- function (comp_base, comp_intv){
  merge_dat <- total_burden [comp == comp_base][total_burden [comp == comp_intv],
                                                on = .(country_name, country)]
  merge_dat <- merge_dat [, .(comp_set = paste0 (i.comp, "_VS_", comp),
                              country_name, country,
                              avt_cases  = total_cases - i.total_cases,
                              avt_deaths = total_deaths - i.total_deaths,
                              avt_dalys  = total_dalys - i.total_dalys,
                              add_doses  = i.total_doses - total_doses)]
  merge_dat [, nnv := add_doses/avt_cases]
  return(merge_dat)
}

all_avtnnv <- rbind (cal_avtnnv ("nomcv", "mcv1"),
                     cal_avtnnv ("mcv1", "mcv1-mcv2"),
                     cal_avtnnv ("mcv1", "mcv1-sia"),
                     cal_avtnnv ("mcv1-sia", "mcv1-mcv2-sia"),
                     cal_avtnnv ("mcv1-mcv2", "mcv1-mcv2-sia"),
                     cal_avtnnv ("mcv1", "mcv1-mcv2alt"),
                     cal_avtnnv ("mcv1-mcv2alt", "mcv1-mcv2alt-sia"),
                     cal_avtnnv ("mcv1-sia", "mcv1-mcv2alt-sia"))
output_case <- setDT (pivot_wider (all_avtnnv [, c("comp_set", "country_name", "country", "avt_cases")],
                                        values_from = avt_cases,
                                        names_from = comp_set))
output_case [, country := factor (country, levels = c(eva_ctries, "Global"))]
setorder (output_case, country)

output_death <- setDT (pivot_wider (all_avtnnv [, c("comp_set", "country_name", "country", "avt_deaths")],
                                    values_from = avt_deaths,
                                    names_from = comp_set))
output_death [, country := factor (country, levels = c(eva_ctries, "Global"))]
setorder (output_death, country)

output_daly <- all_avtnnv [, c("comp_set", "country_name", "country", "avt_dalys")]
output_daly [, avt_dalys_K := avt_dalys/1000]
output_daly <- setDT (pivot_wider (output_daly [, !c("avt_dalys")],
                                   values_from = avt_dalys_K,
                                   names_from = comp_set))
output_daly [, country := factor (country, levels = c(eva_ctries, "Global"))]
setorder (output_daly, country)

output_nnv <- setDT (pivot_wider (all_avtnnv [, c("comp_set", "country_name", "country", "nnv")],
                                  values_from = nnv,
                                  names_from = comp_set))
output_nnv [, country := factor (country, levels = c(eva_ctries, "Global"))]
setorder (output_nnv, country)
# fwrite (x = output_nnv, file = "nnv.csv")

# adjust the order for manuscript
output_tab1_case <- output_case [, 1:6]
names(output_tab1_case)[3:6] <- paste0 (names(output_tab1_case)[3:6], "_case")
output_tab1_death <- output_death [, 3:6]
names(output_tab1_death) <- paste0 (names(output_tab1_death), "_death")
output_tab1 <- cbind (output_tab1_case, output_tab1_death)
output_tab1 <- output_tab1 [, c(1:3,7,4,8,5,9,6,10)]

value_cols <- names(output_tab1)[3:10]
output_tab1 [, (value_cols) := lapply (.SD, function (avt_burden){
                                         return (avt_burden/1000)}),
             .SDcols = value_cols ]
fwrite (x = output_tab1, file = "tab1_avtcasedeath.csv")


# ------------------------------------------------------------------------------
## plot overall cases and incidence trends
# ------------------------------------------------------------------------------
plt_measure <- function (sel_measure, sel_scns, sel_ylab, add.WHOcase){
  plt <- ggplot (data = cum_burden [comp %in% sel_scns],
                 aes(x = year, y = get(sel_measure), colour = comp)) +
    geom_line (size = 0.9, alpha = 0.7) +
    facet_wrap (vars(country_name), scales = "free", ncol = 5) +
    labs (x = "year", y = sel_ylab) +
    scale_colour_brewer("Scenarios", palette = "Dark2") +
    scale_y_log10 () +
    theme_bw () +
    theme (legend.position = "bottom") # legend.position = c(0.93, 0.75)
  if (add.WHOcase) {
      plt <- plt + geom_point (data = data_WHOcase,
                           aes(x = year, y = notifs, shape = "WHO reported cases"),
                           size = 2, colour = "grey50") +
    scale_shape_manual ("", values = 4)
  }
  return(plt)
}
ggsave ("plot/case-overview.pdf",
        plt_measure("cases", vac_stgs[c(1:3,5:6)], "Number of measles cases", T),
        height = 9, width = 15)

ggsave ("plot/case-overview-siaplan.pdf",
        plt_measure("cases", vac_stgs[c(6,8)], "Number of measles cases", T),
        height = 9, width = 15)

ggsave ("plot/incidence-overview.pdf",
        plt_measure("incrateM", vac_stgs[c(1:3,5:6)], "Measles incidence rate (per million)", F),
        height = 9, width = 15)


# ------------------------------------------------------------------------------
## plot incidence trends for a single scenario
# ------------------------------------------------------------------------------
plt_incrate <- function (sel_scn){
  plt <- ggplot (data = cum_burden [comp == sel_scn],
                 aes(x = year, y = incrateM, colour = country_name)) +
    geom_line (size = 0.9, alpha = 0.7) +
    labs (x = "year", y = "Measles incidence rate (per million)") +
    scale_colour_discrete("Country") +
    theme_bw () +
    theme (legend.position = "bottom")
  return(plt)
}
plt_incrate (vac_stgs[6])


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
                                       labels = c("zero-dose", "single-dose", "multi-dose")))]

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
    geom_col (aes(fill = measure), width = 0.8, colour = "white", size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scale = "free_y", ncol = 5) +
    labs (x = "", y = "Number of doses over 2000-2020 (millions)") +
    scale_fill_manual ("Population reached",
                       values = c("#00468b", "#ed0000", "#42b540")) +
    theme_bw () +
    theme (legend.position = "top",
           legend.text = element_text (size = 10),
           axis.text.x = element_text (angle = 60, hjust = 1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0.2, 0.2, 0, 0.4, "cm"))
  return(plt)
}
pdf ("plot/dose-sum.pdf", width = 10, height = 8)
print (plt_dose_sum (eva_ctries, vac_stg_names[c(2,3,5,6)]))
dev.off()


# ------------------------------------------------------------------------------
## plot averted cases and NNV for sensitivity analysis
# ------------------------------------------------------------------------------
all_avtnnv_senanl <- all_avtnnv [comp_set %in% c("mcv1-mcv2_VS_mcv1",
                                                 "mcv1-sia_VS_mcv1",
                                                 "mcv1-mcv2alt_VS_mcv1")]
all_avtnnv_senanl [, `:=` (country_name = factor (country_name,
                                                  levels = country_names[eva_ctries]),
                           comp_set = factor (comp_set,
                                              labels = c("MCV1 + MCV2",
                                                         "MCV1 + MCV2 (alternative)",
                                                         "MCV1 + SIAs")),
                           avt_cases_M = ifelse (avt_cases <= 0.0001, NA, avt_cases/1e6))]


plt_senanl <- function (sel_mea, sel_ylab, add.legend){
  plt <- ggplot (data = all_avtnnv_senanl [country != "Global"],
                 aes(x = country_name, y = get(sel_mea))) +
    geom_col (aes(fill = comp_set), width = 0.8, colour = NA, size = 0,
              position = position_dodge(), colour = "white") +
    labs (x = "", y = sel_ylab) +
    scale_fill_manual ("Scenario compared to MCV1 only",
                       values = c("#00468b","#0099b4", "#ed0000")) +
    theme_bw () +
    theme (legend.position = ifelse (add.legend, "top", "none"),
           legend.text = element_text (size = 9.5),
           axis.text.x = element_text (angle = 60, hjust = 1, size = 10.2),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0.2, 0.2, 0, 0.7, "cm"))
  return(plt)
}

ggsave ("plot/nnv_avtcase_senanl.pdf",
        ggarrange (plt_senanl ("avt_cases_M", "Averted cases (millions)", T),
                   plt_senanl ("nnv", "Number needed to vaccinate", F),
                   common.legend = T, ncol = 1, vjust = -1.2,
                   labels = c("A", "B"), heights = c(5, 4.5)),
        height = 6, width = 11)

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
  subset_data [, `:=` (measure = factor (measure, labels = c("zero-dose", "single-dose", "multi-dose")))]

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
                                     labels = c("zero-dose", "single-dose", "multi-dose"))]

  plt <- ggplot (data = plt_data_sum [measure_type == sel_mea &
                                        comp %in% sel_scns &
                                        value > 0], # remove bar borders in the plot for zero values
                 aes(x = comp, y = value/1e6, group = measure)) +
    geom_col (aes(fill = measure), width = 0.8, colour = NA, size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scales = "free_y", ncol = 5) +
    labs (x = " ", y = sel_ylab) +
    scale_fill_manual ("MCV state of population",
                       values = c("#00468b", "#ed0000", "#42b540")) +
    theme_bw () +
    theme (legend.position = "top",
           legend.text = element_text (size = 10),
           axis.text.x = element_text (angle = 60, hjust = 1),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0.2, 0.2, 0, 0.4, "cm"))
  return(plt)
}

pdf ("plot/case-by-dose-sum-2016to2020.pdf", width = 10, height = 8)
print (plt_burden_byDose_sum (0, 2016:2020, "Cases", vac_stg_names[c(2,3,5,6)],
                              "Oranges", "Cases over 2000-2020 (millions)"))
dev.off()

pdf ("plot/case-by-dose-sum-2016to2020-age2up.pdf", width = 10, height = 8)
print (plt_burden_byDose_sum (2, 2016:2020, "Cases", vac_stg_names[c(2,3,5,6)],
                              "Oranges", "Cases over 2016-2020 (millions), age 2+"))
dev.off()

# pdf ("plot/death-by-dose-sum.pdf", width = 12, height = 12)
# print (plt_burden_byDose_sum (0, "Deaths", vac_stg_names[c(2,3,5,6)], "Purples"))
# dev.off()
