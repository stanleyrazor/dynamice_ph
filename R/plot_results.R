# plot_results.R
# update: 2022/06/14

library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(stringr)
library(readxl)
library(tidyr)

rm (list = ls())
vac_stgs <- c("nomcv",               # (1) no vaccination
              "mcv1",                # (2) MCV1 only
              "mcv1-mcv2",           # (3) MCV1 + MCV2
              "mcv1-mcv2-sia",       # (4) MCV1 + MCV2 + SIA
              "mcv1-sia",            # (5) MCV1 + SIA
              "mcv1-mcv2alt",        # (6) MCV1 + MCV2(early intro)
              "mcv1-mcv2alt-sia",    # (7) MCV1 + MCV2(early intro) + SIA
              "mcv1-mcv2-siaalt1",   # (8) MCV1 + MCV2 + SIA(zero dose first)
              "mcv1-mcv2-siaalt2",   # (9) MCV1 + MCV2 + SIA(already vaccinated first)
              "mcv1-siaalt1",        # (10) MCV1 + SIA(zero dose first)
              "mcv1-siaalt2",        # (11) MCV1 + SIA(already vaccinated first)
              "mcv1-mcv2alt-siaalt1" # (12) MCV1 + MCV2(early intro) + SIA(zero dose first)
)
vac_stg_names <- c("No vaccination",
                   "MCV1",
                   "MCV1 + MCV2",
                   "MCV1 + MCV2 + SIAs",
                   "MCV1 + SIAs",
                   "MCV1 + MCV2 (early intro)",
                   "MCV1 + MCV2 (early intro) + SIAs",
                   "MCV1 + MCV2 + SIAs (zero-dose first)",
                   "MCV1 + MCV2 + SIAs (vaccinated first)",
                   "MCV1 + SIAs (zero-dose first)",
                   "MCV1 + SIAs (vaccinated first)",
                   "MCV1 + MCV2 (early intro) + SIAs (zero-dose first)")

eva_ctries <- c("IND", "NGA", "IDN", "ETH", "CHN",
                "PHL", "UGA", "COD", "PAK", "AGO",
                "MDG", "UKR", "MWI", "SOM")

custom_palette <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                    "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#00468B99",
                    "#ED000099", "#42B54099", "#0099B499", "#925E9F99", "#FDAF9199",
                    "#AD002A99", "#ADB6B699", "#1B191999", "#00468B66", "#ED000066")


# ------------------------------------------------------------------------------
## load and combine model outputs
# ------------------------------------------------------------------------------
# burden and vaccine doses
file_burden  <- NULL
for (scname in vac_stgs){
  scn_burden <- fread (paste0 (getwd(), "/central_burden_estimate/Portnoy/", #"/previous_res/",
                               "central_burden_estimate_", scname, ".csv"))
  scn_burden [, `:=` (cases = cases0d + cases1d + cases2d,
                      deaths = deaths0d + deaths1d + deaths2d)]
  scn_burden [, comp := scname]
  file_burden <- rbind (file_burden, scn_burden)
}
file_burden [country == "COD", country_name := "DRC"]

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
## plotting style functions
# ------------------------------------------------------------------------------
# add a blank window
plt_blank <- ggplot() + geom_blank(aes(1,1)) +
  theme (plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank())

# show non-scientific numbers for axis and allow decimals for small values only
format_num_plain <- function (x, ...) {
  format (x, ..., scientific = FALSE, drop0trailing = TRUE)
}


# ------------------------------------------------------------------------------
## plot country burden over time
# ------------------------------------------------------------------------------
pltdata_ctry_burden <- cum_burden [comp %in% vac_stgs[c(1,2,3,5,4)],
                                   .(comp, year, country_name, cases, deaths, dalys)]
pltdata_ctry_burden <- setDT (pivot_longer (pltdata_ctry_burden,
                                            cols = cases:dalys,
                                            names_to = "measure",
                                            values_to = "value"))
pltdata_ctry_burden [, `:=` (measure = factor (measure, levels = c("cases", "deaths", "dalys"),
                                               labels = c("Cases", "Deaths", "DALYs")),
                             country_name = factor (country_name, levels = country_names[eva_ctries]),
                             comp = factor (comp, levels = vac_stgs[c(1,2,3,5,4)],
                                            labels = vac_stg_names[c(1,2,3,5,4)]))]

pdf("plot/figS3_burden-trend.pdf", height = 8, width = 14)
ggplot (data = pltdata_ctry_burden,
        aes(x = year, y = value/1e6, fill = country_name)) +
  geom_area () +
  facet_grid (rows = vars(measure), cols = vars(comp), scales = "free_y") +
  labs (x = "Year", y = "Estimated health burden (millions)", fill = "Country") +
  scale_fill_manual (values = custom_palette) +
  #scale_x_continuous (labels = label_number (accuracy = 1)) +
  theme_bw () +
  theme (legend.text = element_text (size = 14),
         legend.title = element_text (size = 15),
         axis.title.x = element_text (size = 14, margin = margin (t = 10)),
         axis.title.y = element_text (size = 14, margin = margin (r = 10)),
         axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
         axis.text.y = element_text (size = 10),
         strip.text.x = element_text (size = 14),
         strip.text.y = element_text (size = 14),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         plot.margin = margin (0.2, 0.2, 0.2, 0.2, "cm"))
dev.off()


# ------------------------------------------------------------------------------
## plot incidence trends by scenarios
# ------------------------------------------------------------------------------
# adjust scenario order to allow 'MCV1 only' at the top
pltdat_incrateM <- cum_burden [comp %in% vac_stgs[c(1,2,3,5,4)]]
pltdat_incrateM [, `:=` (country_name = factor (country_name,
                                                levels = country_names[eva_ctries]),
                         comp = factor (comp, levels = vac_stgs [c(4,5,3,2,1)],
                                        labels = vac_stg_names [c(4,5,3,2,1)]))]

pdf (file = "plot/fig2_incrateM.pdf", width = 14, height = 7)
ggplot (data = pltdat_incrateM,
        aes(x = year, y = incrateM, colour = comp)) +
  geom_line (size = 0.9, alpha = 1) +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  labs (x = "Year", y = "Estimated measles incidence rate \nper million population") +
  scale_colour_manual ("Delivery strategies",
                       values = c(custom_palette[c(5,2,1,3)], "grey 50")) +
  # scale_y_log10 (labels = format_num_plain) +
  theme_bw () +
  theme (legend.position = c(0.9, 0.15),
         legend.text = element_text (size = 12),
         legend.title = element_text (size = 12),
         axis.title.x = element_text (size = 15, margin = margin(t = 15)),
         axis.title.y = element_text (size = 15, margin = margin(r = 15)),
         axis.text.x = element_text (size = 10.5),
         axis.text.y = element_text (size = 10.5),
         strip.text.x = element_text (size = 12)) +
  guides (colour = guide_legend (reverse = TRUE))
# geom_point (data = data_WHOcase,
#             aes(x = year, y = incrateM, shape = "WHO reported cases"),
#             size = 2, colour = "gold") +
#   scale_shape_manual ("", values = 4)
dev.off()


# Compare the 'optimal-use' scenarios
pltdat_incrateM_opt <- cum_burden [comp %in% vac_stgs[c(4,7,8,12)]]
pltdat_incrateM_opt [, `:=` (country_name = factor (country_name,
                                                    levels = country_names[eva_ctries]),
                             comp = factor (comp, levels = vac_stgs [c(4,7,8,12)],
                                            labels = vac_stg_names [c(4,7,8,12)]))]

pdf (file = "plot/figS5_incrateM_optim.pdf", width = 14, height = 8)
ggplot (data = pltdat_incrateM_opt,
        aes(x = year, y = incrateM, colour = comp)) +
  # geom_point (data = data_WHOcase,
  #             aes(x = year, y = incrateM, shape = "WHO notifications"),
  #             size = 2, colour = "goldenrod3") +
  geom_line (size = 0.9, alpha = 1) +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  scale_colour_manual ("Delivery strategies",
                       values = custom_palette[c(5,4,6,8)]) +
  scale_shape_manual ("", values = 4) +
  scale_y_log10 (labels = format_num_plain) +
  labs (x = "Year", y = "Estimated measles incidence rate \nper million population (log scale)") +
  theme_bw () +
  theme (legend.position = "bottom",
         legend.text = element_text (size = 11.5),
         legend.title = element_text (size = 12),
         plot.margin = unit (c(0.1, 0.25, 0.1, 0.1), "cm"),
         axis.title.x = element_text (size = 15, margin = margin(t = 15)),
         axis.title.y = element_text (size = 15, margin = margin(r = 10)),
         axis.text.x = element_text (size = 11),
         axis.text.y = element_text (size = 11),
         strip.text.x = element_text (size = 12.5)) +
  guides (colour = guide_legend (order = 1, ncol = 2, byrow = T))
  # guides (colour = guide_legend (order = 1, ncol = 2, byrow = T),
  #         shape = guide_legend (order = 2))
dev.off()


# ------------------------------------------------------------------------------
## plot susceptible population and birth cohort by scenarios
# ------------------------------------------------------------------------------
# birth cohort: not changed by scenario
pltdat_sus_age0 <- copy (file_burden [age == 0 & comp == vac_stgs[1],
                                      c("country_name", "year", "comp", "pops")]) [, comp := "birth"]
pltdat_sus_age0  [, country_name := factor (country_name,
                                            levels = country_names[eva_ctries])]
pltdat_sus <- file_burden [age < 5 & comp %in% vac_stgs[c(1,2,3,5,4)],
                           lapply (.SD, sum), .SDcols = "popsSus",
                           by = c("country", "year", "country_name", "comp")]
pltdat_sus [, `:=` (country_name = factor (country_name,
                                           levels = country_names[eva_ctries]),
                    comp = factor (comp, levels = c(vac_stgs[c(4,5,3,2,1)]),
                                   labels = c(vac_stg_names[c(4,5,3,2,1)])))]
# plot
pdf (file = "plot/fig3_sus-birth.pdf", width = 14, height = 7)
ggplot () +
  geom_line (data = pltdat_sus[comp != "Birth cohort"],
             aes(x = year, y = popsSus/1e6, colour = comp),
             size = 0.9, alpha = 1, linetype = 1) +
  geom_line (data = pltdat_sus_age0,
             aes(x = year, y = pops/1e6, linetype = "birth"),
             size = 0.8, alpha = 0.85) +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  scale_colour_manual ("Delivery strategies",
                       values = c(custom_palette[c(5,2,1,3)], "grey 50")) +
  scale_linetype_manual ("", values = c("birth" = 3), labels = "Annual birth cohort") +
  labs (x = "Year", y = "Estimated number of susceptible population < 5 y/o (millions)") +
  theme_bw () +
  theme (legend.position = c(0.9, 0.1),
         legend.text = element_text (size = 12),
         legend.title = element_text (size = 12),
         axis.title.x = element_text (size = 15, margin = margin(t = 15)),
         axis.title.y = element_text (size = 15, margin = margin(r = 15)),
         axis.text.x = element_text (size = 10.5),
         axis.text.y = element_text (size = 10.5),
         strip.text.x = element_text (size = 12)) +
  guides (colour = guide_legend (reverse = TRUE))
dev.off()

# calculate number of years with susceptibles > birth cohort
compdat_sus <- pltdat_sus [pltdat_sus_age0 [, !c("comp")], on = c("country_name", "year")]
compdat_sus [, outbreak := ifelse(popsSus >= pops, 1, 0)]
tab_sus <- compdat_sus [, .(pr_outbreak = 1-sum(outbreak)/21), by = c("country", "comp")]

tab_sus <- setDT (pivot_wider (tab_sus, values_from = pr_outbreak, names_from = comp))
tab_sus [, country := factor (country, levels = eva_ctries)]
setorder (tab_sus, country)
fwrite (x = tab_sus [, c(1:4,6,5)], file = "tabs3_yr-sus-outbreak.csv")

tab_sus [, lapply(.SD, mean), .SDcols = 3:6]


# ------------------------------------------------------------------------------
## calculate averted burden and number needed to vaccinate (NNV)
# ------------------------------------------------------------------------------
# get averted burden, additional doses, NNV, and proportion of reduction
cal_avtnnv <- function (cum_burden, comp_base, comp_intv){
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

  merge_dat <- total_burden [comp == comp_base][total_burden [comp == comp_intv],
                                                on = .(country_name, country)]
  merge_dat <- merge_dat [, .(comp_set = paste0 (i.comp, "_VS_", comp),
                              country_name, country,
                              avt_cases  = total_cases - i.total_cases,
                              avt_deaths = total_deaths - i.total_deaths,
                              avt_dalys  = total_dalys - i.total_dalys,
                              add_doses  = i.total_doses - total_doses,
                              pr_red_cases = (total_cases - i.total_cases)/total_cases)]
  merge_dat [, nnv := add_doses/avt_cases]
  return(merge_dat)
}

# check absolute and relative case reduction
sel_vacc_impact <- cal_avtnnv (cum_burden, "nomcv", "mcv1")
setorder (sel_vacc_impact, pr_red_cases)
sel_vacc_impact [country != "Global"]
setorder (sel_vacc_impact, avt_cases)
sel_vacc_impact [country != "Global"]
sel_vacc_impact <- cal_avtnnv (cum_burden, "mcv1-mcv2", "mcv1-mcv2alt")
# Global averted cases: 96804681

# combine results of all comparison pairs
all_avtnnv <- rbind (cal_avtnnv (cum_burden, "nomcv", "mcv1"),
                     cal_avtnnv (cum_burden, "nomcv", "mcv1-sia"),
                     cal_avtnnv (cum_burden, "nomcv", "mcv1-mcv2"),
                     cal_avtnnv (cum_burden, "nomcv", "mcv1-mcv2-sia"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-sia"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2"),
                     cal_avtnnv (cum_burden, "mcv1-sia", "mcv1-mcv2-sia"),
                     cal_avtnnv (cum_burden, "mcv1-mcv2", "mcv1-mcv2-sia"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2alt"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-siaalt1"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-siaalt2"))
                     # cal_avtnnv (cum_burden, "mcv1-mcv2alt", "mcv1-mcv2alt-sia"),
                     # cal_avtnnv (cum_burden, "mcv1-sia", "mcv1-mcv2alt-sia"))

# averted cases, deaths, and dalys
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


# Table 1: Averted cases
# country order by MCV2 introduction
eva_ctries_mcv2 <- c("IND", "IDN", "CHN", "PHL", "PAK", "AGO", "UKR", "MWI",
                     "NGA","ETH", "MDG",
                     "UGA","COD", "SOM")
# output_tab1_case <- output_case [, 1:7]
# names(output_tab1_case)[3:7] <- paste0 (names(output_tab1_case)[3:7], "_case")
# output_tab1_death <- output_death [, 3:7]
# names(output_tab1_death) <- paste0 (names(output_tab1_death), "_death")
# output_tab1 <- cbind (output_tab1_case, output_tab1_death)
# output_tab1 <- output_tab1 [, c(1:3,8,4,9,5,10,6,11,7,12)]

value_cols <- names(output_case)[3:8]
output_case [, (value_cols) := lapply (.SD, function (avt_burden){
  return (avt_burden/1000)}),
  .SDcols = value_cols]
output_case [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
setorder (output_case, country)
fwrite (x = output_case[, 1:8], file = "tab1_avtcase.csv")

# Table S2: Averted deaths
output_death [, (value_cols) := lapply (.SD, function (avt_burden){
  return (avt_burden/1000)}),
  .SDcols = value_cols]
output_death [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
setorder (output_death, country)
fwrite (x = output_death[, 1:8], file = "tabs2_avtdeath.csv")

# Table S3: Averted DALYs
output_daly [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
setorder (output_daly, country)
fwrite (x = output_daly[, 1:8], file = "tabs3_avtdaly.csv")

# Table 2: NNV
# not calculated for countries have not yet introduced MCV2 when MCV1 only is the comparator
output_nnv <- setDT (pivot_wider (all_avtnnv [comp_set %in% c("mcv1_VS_nomcv",
                                                              "mcv1-sia_VS_mcv1",
                                                              "mcv1-mcv2_VS_mcv1",
                                                              "mcv1-mcv2-sia_VS_mcv1-mcv2",
                                                              "mcv1-mcv2-sia_VS_mcv1-sia"),
                                              c("comp_set", "country_name", "country", "nnv")],
                                  values_from = nnv,
                                  names_from = comp_set))
output_nnv [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
setorder (output_nnv, country)
fwrite (x = output_nnv, file = "tab2_nnv.csv")

output_nnv [country_name != "Global",
            lapply (.SD, function(x) median (x [x<1000], na.rm = T)), .SDcols = 3:7] # avoid NNV with small-size denominator


# ------------------------------------------------------------------------------
## plot dose distribution
# ------------------------------------------------------------------------------
# process data for doses
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
                     measure = factor (measure,
                                       levels = c("doses0", "doses1", "doses2"),
                                       labels = c("zero-dose", "single-dose", "multi-dose")))]

# # plot dose by time
# plt_dose <- function (sel_ctries){
#   plt <- ggplot (data = pltdat_dose [country %in% sel_ctries & value > 0], # remove bar borders in the plot for zero values
#                  aes(x = year, y = value/1e6)) +
#     geom_col (aes(fill = measure), width = 1, colour = NA, size = 0,
#               position = position_stack (reverse = TRUE)) +
#     facet_grid (rows = vars(country_name), cols = vars(comp), scales = "free_y") +
#     labs (x = "year", y = "Number of doses (millions)") +
#     scale_fill_brewer ("Population reached", palette = "Set2") +
#     theme_bw ()
#   return(plt)
# }
#
# pdf ("plot/dose.pdf", width = 14, height = 6)
# for (ig in 0:3){
#   print (plt_dose (eva_ctries[ig*5 + (1:5)]))
# }
# dev.off()

# plot cumulative doses over 2000-2020 by scenarios
pltdat_dose_sum <- pltdat_dose [, .(total_dose = sum(value)),
                                by = c("comp", "country_name", "country", "measure")]
plt_dose_sum <- function (sel_ctries, sel_scns){
  pltdat_tmp <- pltdat_dose_sum [country %in% sel_ctries & comp %in% vac_stgs[sel_scns]
                                 & total_dose > 0] # remove bar borders in the plot for zero values
  pltdat_tmp [, comp := factor (comp, levels = vac_stgs[sel_scns],
                                labels = vac_stg_names[sel_scns])]
  plt <- ggplot (data = pltdat_tmp, aes(x = comp, y = total_dose/1e6)) +
    geom_col (aes(fill = measure), width = 0.8, colour = "white", size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scale = "free_y", ncol = 5) +
    labs (x = "", y = "Estimated number of doses over 2000-2020 (millions)") +
    scale_fill_manual ("Predicted vaccination state of population reached",
                       values = c("#42b540", "#00468b", "#ed0000")) +
    theme_bw () +
    theme (legend.position  = "top",
           legend.text = element_text (size = 14),
           legend.title = element_text (size = 15),
           axis.title.y = element_text (size = 15, vjust = 2),
           axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           strip.text.x = element_text (size = 12),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0.1, 0.2, 0, 0.6, "cm"))
  return(plt)
}
pdf ("plot/fig4_dose-sum.pdf", width = 11, height = 7)
print (plt_dose_sum (eva_ctries, c(2,3,5,4)))
dev.off()


# ------------------------------------------------------------------------------
## plot averted cases and NNV for sensitivity analysis
# ------------------------------------------------------------------------------
all_avtnnv_senanl <- all_avtnnv [comp_set %in% c("mcv1-mcv2_VS_mcv1",
                                                 "mcv1-mcv2alt_VS_mcv1",
                                                 "mcv1-sia_VS_mcv1",
                                                 "mcv1-siaalt1_VS_mcv1",
                                                 "mcv1-siaalt2_VS_mcv1")]
all_avtnnv_senanl [, `:=` (country_name = factor (country_name,
                                                  levels = country_names[eva_ctries]),
                           comp_set = factor (comp_set,
                                              levels = c("mcv1-mcv2_VS_mcv1",
                                                         "mcv1-mcv2alt_VS_mcv1",
                                                         "mcv1-sia_VS_mcv1",
                                                         "mcv1-siaalt1_VS_mcv1",
                                                         "mcv1-siaalt2_VS_mcv1"),
                                              labels = c("MCV1 + MCV2",
                                                         "MCV1 + MCV2 (early intro)",
                                                         "MCV1 + SIAs",
                                                         "MCV1 + SIAs (zero-dose first)",
                                                         "MCV1 + SIAs (vaccinated first)")),
                           avt_cases_M = ifelse (avt_cases <= 0.0001, NA, avt_cases/1e6))]

# # MCV2 early introduction: calculate proportion of increase in averted cases
# mcv2_avtnnv_senanl <- all_avtnnv_senanl [, c("comp_set", "country","avt_cases")]
# mcv2_avtnnv_senanl <- copy (mcv2_avtnnv_senanl [comp_set == "MCV1 + MCV2"])[mcv2_avtnnv_senanl [comp_set == "MCV1 + MCV2 (early intro)"], on = "country"]
# mcv2_avtnnv_senanl [, pr_inc := (i.avt_cases-avt_cases)/avt_cases]
# setorder(mcv2_avtnnv_senanl, -pr_inc)
# mcv2_avtnnv_senanl [country != "Global"]
#
# median (all_avtnnv_senanl [comp_set == "MCV1 + MCV2 (early intro)" & country != "Global", nnv])
#
# # SIA zero-dose first: calculate proportion of increase in averted cases
# sia1_avtnnv_senanl <- all_avtnnv_senanl [, c("comp_set", "country","avt_cases")]
# sia1_avtnnv_senanl <- copy (sia1_avtnnv_senanl [comp_set ==  "MCV1 + SIAs"])[sia1_avtnnv_senanl [comp_set == "MCV1 + SIAs (zero-dose first)"], on = "country"]
# sia1_avtnnv_senanl [, pr_inc := (i.avt_cases-avt_cases)/avt_cases]
# setorder(sia1_avtnnv_senanl, -pr_inc)
# sia1_avtnnv_senanl [country != "Global"]
#
# median (all_avtnnv_senanl [comp_set == "MCV1 + SIAs (zero-dose first)" & country != "Global", nnv])
# median (all_avtnnv_senanl [comp_set == "MCV1 + SIAs (vaccinated first)" & country != "Global", nnv])
# all_avtnnv_senanl [comp_set == "MCV1 + SIAs (zero-dose first)" & country != "Global", nnv] - all_avtnnv_senanl [comp_set == "MCV1 + MCV2 (early intro)" & country != "Global", nnv]


# plot
plt_senanl <- function (sel_mea, sel_ylab, sel_lgdpos){
  plt <- ggplot (data = all_avtnnv_senanl [country != "Global"],
                 aes(x = country_name, y = get(sel_mea))) +
    geom_col (aes(fill = comp_set), width = 0.8, colour = NA, size = 0,
              position = position_dodge()) +
    labs (x = "", y = sel_ylab) +
    scale_fill_manual ("Delivery strategy compared to MCV1 only",
                       values = custom_palette [c(1,4,2,6,7)]) +
    theme_bw () +
    theme (legend.position = sel_lgdpos,
           legend.text = element_text (size = 13),
           legend.title = element_text (size = 14),
           axis.title.y = element_text (size = 14, vjust = 2),
           axis.text.x = element_text (size = 12, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0, 0.2, 0, 0.6, "cm")) +
    guides (fill = guide_legend (ncol = 3, byrow = T))
  return(plt)
}

ggsave ("plot/fig5_nnv-avtcase-senanl.pdf",
        ggarrange (plt_blank,
                   plt_senanl ("avt_cases_M", "Estimated number of \naverted cases (millions)", c(0.5,1.35)),
                   plt_senanl ("nnv", "Estimated number\nneeded to vaccinate", "none"),
                   ncol = 1, vjust = -1.2,#common.legend = T,
                   labels = c("", "A", "B"), heights = c(1.5, 4, 4)),
        height = 8, width = 11)


# ------------------------------------------------------------------------------
## plot cases and deaths by doses
# ------------------------------------------------------------------------------
# # version 1: by time, scenario, country
# plt_burden_byDose <- function (age_cutoff, sel_yrs, sel_scns, sel_ctries,
#                                sel_title, sel_mea, change_to_pr = F){
#   plt_data <- copy (file_burden [age >= age_cutoff & country %in% sel_ctries &
#                                    year %in% sel_yrs])
#   plt_data <- plt_data [, lapply (.SD, sum), .SDcols = cases0d:deaths2d,
#                         by = c("country", "country_name", "comp", "year")]
#   plt_data [, `:=` (prcases0d = cases0d/(cases0d+cases1d+cases2d),
#                     prcases1d = cases1d/(cases0d+cases1d+cases2d),
#                     prcases2d = cases2d/(cases0d+cases1d+cases2d),
#                     prdeaths0d = deaths0d/(deaths0d+deaths1d+deaths2d),
#                     prdeaths1d = deaths1d/(deaths0d+deaths1d+deaths2d),
#                     prdeaths2d = deaths2d/(deaths0d+deaths1d+deaths2d))]
#   plt_data <- setDT (pivot_longer (plt_data, col = cases0d:prdeaths2d,
#                      names_to = "measure", values_to = "value"))
#   plt_data [, `:=` (country_name = factor (country_name,
#                                            levels = country_names[eva_ctries]),
#                     comp = factor (comp, levels = vac_stgs, labels = vac_stg_names))]
#   subset_data <- plt_data [str_sub(measure,1,5) == sel_mea & comp %in% sel_scns]
#   subset_data [, `:=` (measure = factor (measure, labels = c("zero-dose", "single-dose", "multi-dose")))]
#
#   plt <- ggplot (data = subset_data [value > 0], # remove bar borders in the plot for zero values
#                  aes(x = year, y = value)) +
#     geom_col (aes(fill = measure), width = 1, colour = NA,
#               position = position_stack (reverse = TRUE)) +
#     facet_grid (rows = vars(country_name), cols = vars(comp), scales = "free_y") +
#     labs (x = "year", y = " ", title = sel_title) +
#     scale_fill_brewer ("MCV received", palette = "Oranges") +
#     theme_bw () +
#     theme (legend.position = "right") # legend.position = c(0.93, 0.75)
#
#   if (change_to_pr == T) {
#     plt <- plt + scale_y_continuous (labels = scales::percent)
#   }
#   return(plt)
# }

# pdf ("plot/figS4_case-by-dose.pdf", width = 11, height = 6)
# for (ic in 0:3){
#   print (plt_burden_byDose (0, 2000:2020, vac_stg_names[c(2,5,3,4)],
#                             eva_ctries[ic*5 + (1:5)],
#                             "Estimated number of cases (millions)", "cases", F))
# }
# dev.off()


# version 2: by scenario, country
plt_burden_byDose_sum <- function (age_cutoff, sel_yrs, sel_mea, sel_scns,
                                   sel_palette, sel_ylab){
  plt_data_sum <- copy (file_burden [age >= age_cutoff &
                                       comp %in% vac_stgs[sel_scns] &
                                       year %in% sel_yrs])
  plt_data_sum <- plt_data_sum [, lapply (.SD, sum), .SDcols = cases0d:deaths2d,
                                by = c("country", "country_name", "comp")]
  plt_data_sum <- setDT (pivot_longer (plt_data_sum, col = cases0d:deaths2d,
                                       names_to = "measure", values_to = "value"))
  plt_data_sum [, `:=` (country_name = factor (country_name,
                                               levels = country_names[eva_ctries]),
                        comp = factor (comp, levels = vac_stgs[sel_scns],
                                       labels = vac_stg_names[sel_scns]))]
  plt_data_sum [, measure_type := ifelse (str_sub(measure,1,4) == "case",
                                          "Cases", "Deaths")]
  plt_data_sum [, measure := factor (str_sub (measure,-2,-1),
                                     labels = c("zero-dose", "single-dose", "multi-dose"))]

  plt <- ggplot (data = plt_data_sum [measure_type == sel_mea &
                                        value > 0], # remove bar borders in the plot for zero values
                 aes(x = comp, y = value/1e6, group = measure)) +
    geom_col (aes(fill = measure), width = 0.8, colour = NA, size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scales = "free_y", ncol = 5) +
    labs (x = " ", y = sel_ylab) +
    scale_fill_manual ("Predicted vaccination state of measles cases",
                       values = c("#42b540", "#00468b", "#ed0000")) +
    theme_bw () +
    theme (legend.position  = "top",
           legend.text = element_text (size = 14),
           legend.title = element_text (size = 15),
           axis.title.y = element_text (size = 15, vjust = 2),
           axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           strip.text.x = element_text (size = 12),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0.1, 0.2, 0, 0.6, "cm"))
  return(plt)
}

pdf ("plot/figS4_case-by-dose-sum.pdf", width = 11, height = 7.5)
print (plt_burden_byDose_sum (0, 2000:2020, "Cases", c(2,3,5,4),
                              "Oranges", "Estimated cases over 2000-2020 (millions)"))
dev.off()

# pdf ("plot/case-by-dose-sum-2016to2020-age2up.pdf", width = 10, height = 8)
# print (plt_burden_byDose_sum (2, 2016:2020, "Cases", vac_stg_names[c(2,3,5,6)],
#                               "Oranges", "Cases over 2016-2020 (millions), age 2+"))
# dev.off()


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
## check age distribution of cases and deaths
# ------------------------------------------------------------------------------
global_burden_u5 <- copy (file_burden[age <= 5]) [, lapply (.SD, sum),
                                                  .SDcols = pops:deaths,
                                                  by = "comp"][, agegrp := "under5"]
global_burden_o5 <- copy (file_burden[age > 5]) [, lapply (.SD, sum),
                                                 .SDcols = pops:deaths,
                                                 by = "comp"][, agegrp := "over5"]
global_burden_by_age <- global_burden_u5 [, .(comp, cases, deaths)][global_burden_o5 [,.(comp, cases, deaths)],
                                                                    on = .(comp)]
setnames (global_burden_by_age,
          old = c("cases", "deaths", "i.cases", "i.deaths"),
          new = c("cases_u5", "deaths_u5", "cases_o5", "deaths_o5"))

global_avert_by_age <- global_burden_by_age [comp %in% c("mcv1-mcv2", "mcv1-sia")]
global_avert_by_age [, `:=` (base_cases_u5 = global_burden_by_age [comp == "mcv1", cases_u5],
                             base_cases_o5 = global_burden_by_age [comp == "mcv1", cases_o5],
                             base_deaths_u5 = global_burden_by_age [comp == "mcv1", deaths_u5],
                             base_deaths_o5 = global_burden_by_age [comp == "mcv1", deaths_o5])]
global_avert_by_age [, `:=` (avt_cases_u5 = base_cases_u5 - cases_u5,
                             avt_cases_o5 = base_cases_o5 - cases_o5,
                             avt_deaths_u5 = base_deaths_u5 - deaths_u5,
                             avt_deaths_o5 = base_deaths_o5 - deaths_o5)]
global_avert_by_age [, `:=` (avt_cases_all = avt_cases_u5 + avt_cases_o5,
                             avt_deaths_all = avt_deaths_u5 + avt_deaths_o5)]
global_avert_by_age$avt_cases_u5/global_avert_by_age$avt_cases_all


# ------------------------------------------------------------------------------
## plot relationship between susceptibles and zero-dose population
# ------------------------------------------------------------------------------
pltdat_sus <- file_burden [age < 5, .(prSus = sum(popsSus)/sum(pops),
                                      pr0d  = sum(pops0d)/sum(pops)),
                           by = c("country", "country_name", "year", "comp")]
pltdat_sus [, country_name := factor (country_name, levels = country_names)]

ggplot (data = pltdat_sus [comp != "nomcv"],
        aes(x = pr0d, y = prSus, colour = comp)) +
  geom_point () +
  labs (x = "Proportion of zero-dose", y = "Proportion of susceptibles") +
  scale_colour_brewer ("Scenario", palette = "Set1") +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  theme_bw ()

pdf ("plot/under5sus.pdf", height = 7, width = 11)
ggplot (data = pltdat_sus,
        aes(x = year, y = prSus, colour = comp)) +
  geom_line (size = 0.9, alpha = 0.7) +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  labs (x = "Year", y = "Proportion of susceptible under 5 years old") +
  scale_colour_discrete ("Scenarios") +
  theme_bw () +
  theme (legend.position = "bottom")
dev.off()
