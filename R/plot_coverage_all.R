# plot_coverage_all.R
# generate coverage input files
# follow the file created in process_coverage_sia.R
# update: 2023/04/07

library(data.table)
library(ggplot2)
library(scales)

rm (list = c())
eva_ctries <- c("IND", "NGA", "IDN", "ETH", "CHN",
                "PHL", "UGA", "COD", "PAK", "AGO",
                "MDG", "UKR", "MWI", "SOM")
country_names <- c("India", "Nigeria", "Indonesia", "Ethiopia", "China",
                   "Philippines", "Uganda", "DRC", "Pakistan", "Angola",
                   "Madagascar", "Ukraine", "Malawi", "Somalia")
names(country_names) <- eva_ctries

# folder for saving figures
res_folder <- paste0 (getwd(), "/previous_res/20230401/")


# ------------------------------------------------------------------------------
## plot coverage trends for both routine immunisation and campaigns
# ------------------------------------------------------------------------------
outfile_mcv1_mcv2_sia  <- fread ("coverage/coverage_mcv1-mcv2-sia.csv")
outfile_mcv1_mcv2alt1  <- fread ("coverage/coverage_mcv1-mcv2alt1.csv")
outfile_mcv1_mcv2alt2  <- fread ("coverage/coverage_mcv1-mcv2alt2.csv")
plt_data <- rbind (outfile_mcv1_mcv2_sia,
                   copy (outfile_mcv1_mcv2alt1 [vaccine == "MCV2"])[, vaccine := "MCV2 \n(early intro, fast rollout)"],
                   copy (outfile_mcv1_mcv2alt2 [vaccine == "MCV2"])[, vaccine := "MCV2 \n(early intro, gradual rollout)"])

# update country names
plt_data [country_code == "COD", country := "DRC"]
plt_data [vaccine == "SIA", vaccine := "SIAs"]

# rank countries by IHME burden
plt_data [, country := factor (country, levels = country_names[eva_ctries])]

# figure 1: MCV1, MCV2, SIAs coverage
pdf (paste0 (res_folder, "figures/fig1_covall.pdf"), width = 14, height = 7)
plt_cov <- ggplot (data = plt_data [vaccine %in% c("MCV1", "MCV2") & year >= 2000],
                   aes (x = year, y = coverage, colour = vaccine)) +
  scale_x_continuous (breaks = pretty_breaks ()) +
  geom_line (size = 0.9) +
  facet_wrap (vars(country), ncol = 5) +
  labs (title = " ", x = "Year", y = "Vaccine coverage") +
  theme_bw() +
  theme (legend.position  = c(0.9, 0.08),#"bottom"
         legend.direction = "vertical",
         legend.key.size = unit (1, 'cm'),
         legend.text = element_text (size = 13),
         legend.title = element_text (size = 14),
         axis.title.x = element_text (size = 15, vjust = -0.75),
         axis.title.y = element_text (size = 15, margin = margin(r = 15)),
         axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
         axis.text.y = element_text (size = 10),
         strip.text.x = element_text (size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         plot.margin = margin (0, 0.2, 0.2, 0.2, "cm")) +
  geom_point (data = plt_data [vaccine == "SIAs" & year >= 2000],
              aes (x = year, y = coverage), size = 1.2) +
  scale_colour_manual ("Delivery strategy",
                       values = c("#42b540", "#00468b",  "#ed0000")) +
  guides(color = guide_legend (override.aes = list (shape = c(NA,NA,16),
                                                    linetype = c(1,1,NA),
                                                    size = c(1,1,1.2))))
print(plt_cov)
dev.off()


# figure S1: MCV2, MCV2 (early intro) coverage
pdf (paste0 (res_folder, "figures/figS1_mcv2alt-cov.pdf"), width = 12, height = 7)
ggplot (data = plt_data [vaccine != "SIAs" & year >= 2000],
                        aes (x = year, y = coverage, colour = vaccine, linetype = vaccine)) +
  scale_x_continuous (breaks = pretty_breaks ()) +
  geom_line (size = 1) +
  facet_wrap (vars(country), ncol = 5) +
  labs (title = " ", x = "Year", y = "Vaccine coverage") +
  scale_colour_manual (name = "Delivery strategy",
                       values = c("#42b540", "#00468b", "#0099b4", "#00468B99")) +
  scale_linetype_manual ("Delivery strategy", values = c(1,1,1,2)) +
  theme_bw() +
  theme (legend.position  ="bottom",
         legend.key.size = unit (2, 'cm'),
         legend.text = element_text (size = 13),
         legend.title = element_text (size = 14),
         axis.title.x = element_text (size = 15, vjust = -0.75),
         axis.title.y = element_text (size = 15, margin = margin(r = 15)),
         axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
         axis.text.y = element_text (size = 10),
         strip.text.x = element_text (size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         plot.margin = margin (0, 0.2, 0.2, 0.2, "cm"))
dev.off()

