## ---------------------------
## Purpose: Summarize and plot MEC isolate (EUCAST) MIC and SMG data by species
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

# Load packages
library(patchwork)
library(gtools)
library(data.table)

source("MIC_data_summary.R")

species_colors <- c(
  "#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
  "#882255", "#BBBBBB", "#AA4499", "#DDCC77", "black"
)

species_colors <- species_colors %>%
  set_names(species_count$genus_species)

drugs <- as_labeller(c(
  fluconazole = "Fluconazole",
  micafungin = "Micafungin",
  `amphotericin B` = "Amphotericin B"
))

full_mic <- data.table::copy(mic_info)

# Bin Etest reults with broth MIC levels
full_mic <- full_mic %>%
  mutate(mic50 = replace(mic50, mic50 == ">32", "64")) %>%
  mutate(mic50 = replace(mic50, mic50 == "160", "256"))

full_mic <- full_mic %>%
  mutate(mic90 = replace(mic90, mic90 == "0.023", "0.032"))

full_mic <- full_mic %>%
  mutate(mic50 = case_when((drug == "micafungin" & genus_species == "C. parapsilosis" & mic50 == "2") ~ ">1",
    .default = mic50
  ))

# Drug-specific levels for plotting
flc_levels <- c("0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256")
mcf_amb_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1")

# Drug-specific DFs
mic_flc <- full_mic %>%
  filter(drug == "fluconazole")
mic_flc$mic50 <- factor(mic_flc$mic50, levels = flc_levels)

mic_mcf <- full_mic %>%
  filter(drug == "micafungin")
mic_mcf$mic50 <- factor(mic_mcf$mic50, levels = mcf_amb_levels)

mic_amb <- full_mic %>%
  filter(drug == "amphotericin B")
mic_amb$mic90 <- factor(mic_amb$mic90, levels = mcf_amb_levels)

# Find the series with > 2-fold changes
flc_step_diffs <- mic_flc %>%
  filter(!is.na(series_id)) %>%
  group_by(series_id) %>%
  summarize(diff = max(as.numeric(mic50), na.rm = T) - min(as.numeric(mic50), na.rm = T)) %>%
  filter(diff > 1)

mcf_step_diffs <- mic_mcf %>%
  filter(!is.na(series_id)) %>%
  group_by(series_id) %>%
  summarize(diff = max(as.numeric(mic50), na.rm = T) - min(as.numeric(mic50), na.rm = T)) %>%
  filter(diff > 1)

amb_step_diffs <- mic_amb %>%
  filter(!is.na(series_id)) %>%
  group_by(series_id) %>%
  summarize(diff = max(as.numeric(mic90), na.rm = T) - min(as.numeric(mic90), na.rm = T)) %>%
  filter(diff > 1)

# Plot series AH and BB
bb_flc <- mic_flc %>%
  filter(series_id == "BB")
bb_flc$sample_id <- c("1.1", "1.2", "1.4", "1.3", "2.4", "2.1", "2.2", "2.3")

bb_flc_plot <- bb_flc %>%
  ggplot(aes(x = sample_id, y = mic50, fill = genus_species)) +
  geom_col() +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) +
  theme(axis.text = element_text(size = 8, color = "black")) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  ylab(expression(paste("Fluconazole MIC50, ", mu, "g/mL"))) +
  xlab(NULL)

bb_amb <- mic_amb %>%
  filter(series_id == "BB")
bb_amb$sample_id <- c("1.1", "1.2", "1.4", "1.3", "2.4", "2.1", "2.2", "2.3")

bb_amb_plot <- bb_amb %>%
  ggplot(aes(x = sample_id, y = mic90, fill = genus_species)) +
  geom_col() +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) +
  theme(axis.text = element_text(size = 8, color = "black")) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  ylab(expression(paste("Amphotericin B MIC90, ", mu, "g/mL"))) +
  xlab(NULL)

ah_mcf <- mic_mcf %>%
  filter(series_id == "AH")
ah_mcf$sample_id <- c("1.1", "1.2", "1.3", "1.4", "2.1", "2.2", "3.2", "4.1", "3.1", "4.2", "5.1", "5.2", "6.1", "6.2", "7.1", "7.2", "8.1", "8.2", "9.1", "11.1")

level_order <- mixedsort(ah_mcf$sample_id)

ah_mcf_plot <- ah_mcf %>%
  ggplot(aes(x = sample_id, y = mic50, fill = genus_species)) +
  geom_col(aes(x = factor(sample_id, level = level_order))) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 6, color = "black")) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  ylab(expression(paste("Micafungin MIC50, ", mu, "g/mL"))) +
  xlab(NULL)

ah_amb <- mic_amb %>%
  filter(series_id == "AH")
ah_amb$sample_id <- ah_mcf$sample_id

ah_amb_plot <- ah_amb %>%
  ggplot(aes(x = sample_id, y = mic90, fill = genus_species)) +
  geom_col(aes(x = factor(sample_id, level = level_order))) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 6, color = "black")) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  ylab(expression(paste("Amphotericin B MIC90, ", mu, "g/mL"))) +
  xlab(NULL)

bb_combined <- bb_flc_plot + bb_amb_plot
ah_combined <- ah_mcf_plot + ah_amb_plot

combined_series_diffs <- bb_combined / ah_combined + plot_annotation(tag_levels = list(c("A", "", "B", "")))

ggsave("MEC_serial_MIC_differences.tiff", combined_series_diffs, width = 8, height = 6, units = "in", dpi = 320)
