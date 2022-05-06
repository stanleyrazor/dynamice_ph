# plot_results.R
# update: 2022/05/03

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
              "mcv1-mcv2-sia",     # (4) MCV1 + MCV2 + SIA
              "mcv1-sia",          # (5) MCV1 + SIA
              "mcv1-mcv2alt",      # (6) MCV1 + MCV2(early intro)
              "mcv1-mcv2alt-sia",  # (7) MCV1 + MCV2(early intro) + SIA
              "mcv1-mcv2-siaalt1", # (8) MCV1 + MCV2 + SIA(zero dose first)
              "mcv1-mcv2-siaalt2", # (9) MCV1 + MCV2 + SIA(already vaccinated first)
              "mcv1-siaalt1",      # (10) MCV1 + SIA(zero dose first)
              "mcv1-siaalt2"       # (11) MCV1 + SIA(already vaccinated first)
)
vac_stg_names = c("No vaccination",
                  "MCV1",
                  "MCV1 + MCV2",
                  "MCV1 + MCV2 + SIAs",
                  "MCV1 + SIAs",
                  "MCV1 + MCV2 (early intro)",
                  "MCV1 + MCV2 (early intro) + SIAs",
                  "MCV1 + MCV2 + SIAs (zero-dose first)",
                  "MCV1 + MCV2 + SIAs (vaccinated first)",
                  "MCV1 + SIAs (zero-dose first)",
                  "MCV1 + SIAs (vaccinated first)")

eva_ctries = c("IND", "NGA", "IDN", "ETH", "CHN",
               "PHL", "UGA", "COD", "PAK", "AGO",
               "MDG", "UKR", "MWI", "SOM")


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
custom_palette <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                    "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#00468B99",
                    "#ED000099", "#42B54099", "#0099B499", "#925E9F99", "#FDAF9199",
                    "#AD002A99", "#ADB6B699", "#1B191999", "#00468B66", "#ED000066")

pdf("plot/fig2_burden-trend.pdf", height = 8, width = 14)
ggplot (data = pltdata_ctry_burden,
        aes(x = year, y = value/1e6, fill = country_name)) +
  geom_area () +
  facet_grid (rows = vars(measure), cols = vars(comp), scales = "free_y") +
  labs (x = "Year", y = "Health burden (millions)", fill = "Country") +
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

# check each comparison pair
all_avtnnv <- rbind (cal_avtnnv (cum_burden, "nomcv", "mcv1"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-sia"),
                     cal_avtnnv (cum_burden, "mcv1-sia", "mcv1-mcv2-sia"),
                     cal_avtnnv (cum_burden, "mcv1-mcv2", "mcv1-mcv2-sia"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2alt"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-siaalt1"),
                     cal_avtnnv (cum_burden, "mcv1", "mcv1-siaalt2"),
                     cal_avtnnv (cum_burden, "mcv1-mcv2alt", "mcv1-mcv2alt-sia"),
                     cal_avtnnv (cum_burden, "mcv1-sia", "mcv1-mcv2alt-sia"))

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


# Table 1: Averted cases and deaths
output_tab1_case <- output_case [, 1:7]
names(output_tab1_case)[3:7] <- paste0 (names(output_tab1_case)[3:7], "_case")
output_tab1_death <- output_death [, 3:7]
names(output_tab1_death) <- paste0 (names(output_tab1_death), "_death")
output_tab1 <- cbind (output_tab1_case, output_tab1_death)
output_tab1 <- output_tab1 [, c(1:3,8,4,9,5,10,6,11,7,12)]

value_cols <- names(output_tab1)[3:12]
output_tab1 [, (value_cols) := lapply (.SD, function (avt_burden){
                                         return (avt_burden/1000)}),
             .SDcols = value_cols ]
fwrite (x = output_tab1, file = "tab1_avtcasedeath.csv")

# Table S2: Averted DALYs
fwrite (x = output_daly [, 1:7], file = "tabs2_avtdaly.csv")

# Table 2: NNV
output_nnv <- setDT (pivot_wider (all_avtnnv [, c("comp_set", "country_name", "country", "nnv")],
                                  values_from = nnv,
                                  names_from = comp_set))
output_nnv [, country := factor (country, levels = c(eva_ctries, "Global"))]
setorder (output_nnv, country)
fwrite (x = output_nnv[, c(1,3:7)], file = "tab2_nnv.csv")

output_nnv [country_name != "Global",
            lapply (.SD, function(x) median (x [x<1000], na.rm = T)), .SDcols = 3:7] # avoid NNV with small denominator


# ------------------------------------------------------------------------------
# analyse effect of time periods on NNV
# ------------------------------------------------------------------------------
scn_evaluate <- c("mcv1",  "mcv1-mcv2", "mcv1-sia", "mcv1-mcv2-sia", "mcv1-mcv2-sia")
scn_baseline <- c("nomcv", "mcv1",      "mcv1",     "mcv1-mcv2",     "mcv1-sia")
scn_pairs    <- paste0 (toupper(scn_evaluate), " vs ", toupper(scn_baseline))
anl_timeframe <- list (2000:2020, 2001:2010, 2011:2020)
all_nnv <- NULL
for (it in 1:3){
  total_burden <- file_burden [year %in% anl_timeframe[[it]], lapply (.SD, sum),
                               .SDcols = c("cases", "doses"),
                               by = c("country_name", "comp")]
  for (ip in 1:5){
    merge_dat <- total_burden [comp == scn_baseline[ip]][total_burden [comp == scn_evaluate[ip]],
                                                         on = .(country_name)]
    merge_dat <- merge_dat [, .(comp_set = scn_pairs[ip],
                                country_name,
                                avt_cases  = cases - i.cases,
                                add_doses  = i.doses - doses)]
    merge_dat [, `:=` (nnv = add_doses/avt_cases,
                       time_frame = paste0 (min(anl_timeframe[[it]]), "-", max(anl_timeframe[[it]])))]
    all_nnv <- rbind (all_nnv, merge_dat)
    }
}
all_nnv [, `:=` (country_name = factor (country_name, levels = country_names[eva_ctries]),
                 comp_set = factor (comp_set, levels = scn_pairs))]
all_nnv [avt_cases < 0.5, nnv := NA]
pdf (file = "plot/time-period.pdf", width = 12, height = 9)
ggplot (data = copy(all_nnv)[comp_set %in% scn_pairs[c(1,2,3)]], # [, nnv := ifelse(nnv>50, NA, nnv)], # do not plot nnv > 50
        aes (x = comp_set, y = nnv, colour = time_frame,
             group = paste0 (time_frame, country_name))) +
  geom_point (aes(shape = time_frame), size = 2) +
  geom_line () +
  facet_wrap (vars(country_name), ncol = 5, scales = "free_y") +
  # facet_grid (rows = vars(country_name), cols = vars(comp_set), scales = "free") +
  labs (x = " ", y = "Number needed to vaccinate") +
  scale_colour_discrete ("Time frame") +
  scale_shape (guide = "none") +
  theme_bw () +
  theme (legend.text = element_text (size = 13),
         legend.title = element_text (size = 14),
         axis.title.x = element_text (size = 14, margin = margin (t = 10)),
         axis.title.y = element_text (size = 14, margin = margin (r = 10)),
         axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
         axis.text.y = element_text (size = 10),
         strip.text.x = element_text (size = 12),
         strip.text.y = element_text (size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         plot.margin = margin (0.2, 0.2, 0.2, 0.2, "cm"))
dev.off()


# # ------------------------------------------------------------------------------
# # analyse total effectiveness and administered doses
# # ------------------------------------------------------------------------------
# avert_nomcv <- cum_burden [comp %in% vac_stgs[2:7], lapply (.SD, sum),
#                            .SDcols = c("cases", "doses"),
#                            by = c("country_name", "comp")]
# burden_nomcv <- cum_burden [comp == "nomcv", .(nomcv_case = sum(cases)),
#                             by = c("country_name")]
# avert_nomcv <- avert_nomcv [burden_nomcv, on = .(country_name)]
# avert_nomcv [, `:=` (avt_cases = nomcv_case - cases,
#                      stg_group = ifelse (comp == "mcv1-sia", 2,
#                                          ifelse (comp %in% c("mcv1-mcv2alt", "mcv1-mcv2alt-sia"), 3, 1)))]
# avert_nomcv <- rbind (avert_nomcv,
#                       copy(avert_nomcv [comp %in% c("mcv1", "mcv1-mcv2-sia")])[, stg_group := 2],
#                       copy(avert_nomcv [comp == "mcv1"])[, stg_group := 3])
# avert_nomcv [, `:=` (comp = factor (comp, levels = vac_stgs[c(2,3,5,6,4,7)],
#                                     labels = vac_stg_names[c(2,3,5,6,4,7)]),
#                      country_name = factor (country_name, levels = country_names[eva_ctries]))]
# pdf (file = "plot/average-avtcase-dose.pdf", width = 14, height = 8)
# ggplot (data = avert_nomcv,
#         aes (x = avt_cases/1e6, y = doses/1e6,
#              group = stg_group, colour = as.factor(stg_group))) +
#   geom_point (aes(shape = comp), size = 2.5) +
#   geom_line (size = 0.75, alpha = 0.8) +
#   facet_wrap (vars(country_name), ncol = 5, scales = "free") +
#   labs (x = "Averted cases compared to no vaccination (millions)",
#         y = "Number of doses administered (millions)") +
#   scale_shape_manual ("Delivery strategies", values = 0:5) +
#   scale_colour_discrete ("Order of additional strategies",
#                          labels = c("MCV2 -> SIAs", "SIAs -> MCV2",
#                                     "MCV2 -> SIAs (Alternative)")) +
#   scale_x_continuous (limits = c(0, NA)) +
#   scale_y_continuous (limits = c(0, NA)) +
#   theme_bw () +
#   theme (legend.text = element_text (size = 13),
#          legend.title = element_text (size = 14),
#          axis.title.x = element_text (size = 14, margin = margin (t = 10)),
#          axis.title.y = element_text (size = 14, margin = margin (r = 10)),
#          axis.text.y = element_text (size = 10),
#          strip.text.x = element_text (size = 13),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.background = element_blank(),
#          plot.margin = margin (0.2, 0.2, 0.2, 0.2, "cm"))
# dev.off()
#

# ------------------------------------------------------------------------------
## plot overall cases and incidence trends
# ------------------------------------------------------------------------------
plt_measure <- function (sel_measure, sel_scns, sel_ylab, add.WHOcase){
  plt <- ggplot (data = cum_burden [comp %in% sel_scns],
                 aes(x = year, y = get(sel_measure), colour = comp)) +
    geom_line (size = 0.9, alpha = 0.7) +
    facet_wrap (vars(country_name), scales = "free", ncol = 5) +
    labs (x = "year", y = sel_ylab) +
    # scale_colour_brewer("Scenarios", palette = "Dark2") +
    scale_colour_discrete ("Scenarios") +
    # scale_y_log10 () +
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

ggsave ("plot/incidence-overview.pdf",
        plt_measure("incrateM", vac_stgs[c(1:3,5:6)], "Measles incidence rate (per million)", F),
        height = 8, width = 15)


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
                     comp = factor (comp, levels = vac_stgs, labels = vac_stg_names),
                     measure = factor (measure,
                                       levels = c("doses0", "doses1", "doses2"),
                                       labels = c("zero-dose", "single-dose", "multi-dose")))]

# plot dose by time
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

# plot cumulative doses over 2000-2020 by scenarios
pltdat_dose_sum <- pltdat_dose [, .(total_dose = sum(value)),
                                by = c("comp", "country_name", "country", "measure")]
plt_dose_sum <- function (sel_ctries, sel_scns){
  pltdat_tmp <- pltdat_dose_sum [country %in% sel_ctries & comp %in% sel_scns
                                 & total_dose > 0] # remove bar borders in the plot for zero values
  pltdat_tmp [, comp := factor (comp, levels = sel_scns)]
  plt <- ggplot (data = pltdat_tmp, aes(x = comp, y = total_dose/1e6)) +
    geom_col (aes(fill = measure), width = 0.8, colour = "white", size = 0,
              position = position_stack (reverse = TRUE)) +
    facet_wrap (vars(country_name), scale = "free_y", ncol = 5) +
    labs (x = "", y = "Number of doses over 2000-2020 (millions)") +
    scale_fill_manual ("Vaccination state of population reached",
                       values = c("#00468b", "#ed0000", "#42b540")) +
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
pdf ("plot/fig3_dose-sum.pdf", width = 11, height = 8)
print (plt_dose_sum (eva_ctries, vac_stg_names[c(2,3,5,4)]))
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


plt_senanl <- function (sel_mea, sel_ylab, sel_lgdpos){
  plt <- ggplot (data = all_avtnnv_senanl [country != "Global"],
                 aes(x = country_name, y = get(sel_mea))) +
    geom_col (aes(fill = comp_set), width = 0.8, colour = NA, size = 0,
              position = position_dodge()) +
    labs (x = "", y = sel_ylab) +
    scale_fill_manual ("Scenario compared to MCV1 only",
                       values = custom_palette [c(1,4,2,3,5)]) +
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

ggsave ("plot/fig4_nnv-avtcase-senanl.pdf",
        ggarrange (plt_blank,
                   plt_senanl ("avt_cases_M", "Averted cases\n (millions)", c(0.5,1.35)),
                   plt_senanl ("nnv", "Number needed\n to vaccinate", "none"),
                   ncol = 1, vjust = -1.2,#common.legend = T,
                   labels = c("", "A", "B"), heights = c(1.5, 4, 4)),
        height = 8, width = 11)


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
    scale_fill_manual ("Vaccination state of measles cases",
                       values = c("#00468b", "#ed0000", "#42b540")) +
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

pdf ("plot/figS1_case-by-dose-sum.pdf", width = 11, height = 8)
print (plt_burden_byDose_sum (0, 2000:2020, "Cases", vac_stg_names[c(2,3,5,6)],
                              "Oranges", "Cases over 2000-2020 (millions)"))
dev.off()

pdf ("plot/case-by-dose-sum-2016to2020-age2up.pdf", width = 10, height = 8)
print (plt_burden_byDose_sum (2, 2016:2020, "Cases", vac_stg_names[c(2,3,5,6)],
                              "Oranges", "Cases over 2016-2020 (millions), age 2+"))
dev.off()

# pdf ("plot/death-by-dose-sum.pdf", width = 12, height = 12)
# print (plt_burden_byDose_sum (0, "Deaths", vac_stg_names[c(2,3,5,6)], "Purples"))
# dev.off()
