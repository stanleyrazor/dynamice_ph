# plot_results.R
# update: 2022/01/20

rm (list = c())
library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(stringr)
library(readxl)

vaccine_strategies <- c("nomcv",         # (1) no vaccination
                        "mcv1",          # (2) MCV1 only
                        "mcv1-mcv2",     # (3) MCV1 + MCV2
                        "mcv1-mcv2-sia", # (4) MCV1 + MCV2 + SIA
                        "mcv1-sia",      # (5) MCV1 + SIA
                        "mcv1-mcv2alt"   # (6) MCV1 + MCV2(alternative)
)

# ------------------------------------------------------------------------------
## load and combine model outputs
# ------------------------------------------------------------------------------
# burden and vaccine doses
file_case  <- NULL
for (scname in vaccine_strategies){
  scn_case <- fread( paste0(getwd(), "/central_burden_estimate/Portnoy/",
                            "central_burden_estimate_", scname, ".csv"))
  scn_case [, `:=` (cases = cases0d + cases1d + cases2d,
                    deaths = deaths0d + deaths1d + deaths2d)]
  scn_case [, comp := scname]
  file_case <- rbind(file_case, scn_case)
}

cum_burden <- file_case [, lapply (.SD, sum),
                         .SDcols = cohort_size:deaths,
                         by = c("country", "year", "country_name", "comp")]

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
plt_measure <- function (sel_measure){
  plt <- ggplot (data = cum_burden,
                 aes(x = year, y = get(sel_measure), colour = comp)) +
  geom_line (size = 0.9, alpha = 0.7) +
  facet_wrap (vars(country_name), scales = "free", ncol = 5) +
  labs (title = "Number of measles cases", x = "year", y = " ") +
  geom_point (data = data_WHOcase,
              aes(x = year, y = notifs, shape = "WHO reported cases"),
              size = 2, colour = "grey50") +
  scale_shape_manual ("", values = 4) +
  scale_colour_brewer("Scenarios", palette = "Dark2") +
  theme_bw () +
  theme (legend.position = "bottom") # legend.position = c(0.93, 0.75)
  return(plt)
}
ggsave ("plot/case_overview.pdf", plt_measure("cases"), height = 9, width = 15)

