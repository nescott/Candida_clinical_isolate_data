## ---------------------------
## Purpose: Take MIC_calcs.R values as input, plot Selmecki-style MIC heatmap
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)
## ---------------------------
## load packages
library(tidyverse)
library(patchwork)
library(gtools)
## ---------------------------
# Create empty variable for optional use in MIC_SMG_plots_redcap.R
exclude <- c("")

# Proper sorting of strain IDs for plotting
strain_order <- function(mic_cut) {
  mic_cut %>%
    filter(!strain %in% exclude) -> mic_cut
  mixedsort(mic_cut$strain)
}

# Basic qc plot for data exploration
mic_boxplot <- function(tidy_mic) {
  tidy_mic$strain <- toupper(tidy_mic$strain)
  tidy_mic %>%
    mutate(across(strain, as_factor))
  tidy_mic$strain <- fct_relevel(tidy_mic$strain, mixedsort)
  tidy_mic <- tidy_mic %>%
    filter(concentration != "blank", strain != "blank")
  # tidy_mic$concentration <- fct_rev(tidy_mic$concentration)
  tidy_mic %>%
    ggplot(aes(x = concentration, y = OD600)) +
    facet_wrap(~strain, ncol = 3) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90)) +
    theme_bw() +
    xlab(paste(drug$drug[1], "concentration")) +
    ylab("OD")
}

plotting_coords <- function(final_od, mic_cut, mic_cutpoint, strain_order) {
  final_od %>%
    filter(!strain %in% exclude) -> final_od

  mic_cut %>%
    filter(!strain %in% exclude) %>%
    arrange(factor(strain, levels = strain_order)) -> mic_cut

  mic_cut_x <- setNames(as.numeric(ifelse(mic_cut$mean_norm_OD <= mic_cutpoint, mic_cut$concentration, as.numeric(mic_cut$concentration) + 1)), mic_cut$strain)
  mic_cut_y <- setNames((seq(0.5, n_distinct(final_od$strain), by = 1)), mic_cut$strain)
  mic_cut_yend <- setNames((seq(1.5, (n_distinct(final_od$strain) + 1), by = 1)), mic_cut$strain)
  plotting_coords <- list(
    mic_cut_y = mic_cut_y,
    mic_cut_yend = mic_cut_yend,
    mic_cut_x = mic_cut_x,
    mic_cut_xend = (mic_cut_x)
  )
}

# Plot with strains on x-axis, concentration on y-axis, no labels
mic_plot <- function(final_od, strain_order, plotting_coords) {
  max(final_od$mean_norm_OD) -> scale_lim
  print(scale_lim)
  final_od %>%
    filter(!strain %in% exclude) %>%
    mutate(mean_norm_OD = case_when(mean_norm_OD < 0 ~ 0,
      .default = mean_norm_OD
    )) %>%
    ggplot(aes(y = strain, x = concentration, fill = mean_norm_OD)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(family = "Helvetica", color = "black", size = 10),
      axis.text.y = element_text(family = "Helvetica", color = "black", size = 11, vjust = 0.4),
      legend.text = element_text(family = "Helvetica", color = "black"),
      plot.title = element_text(family = "Helvetica", color = "black", size = 13)
    ) +
    scale_fill_continuous(limits = c(0, scale_lim)) +
    guides(fill = guide_legend(title = "Relative \ngrowth")) +
    ylab(NULL) +
    xlab(NULL) +
    geom_raster(hjust = 1.0) + # adjust tick marks to edge of well
    annotate("segment",
      y = plotting_coords$mic_cut_y,
      yend = plotting_coords$mic_cut_yend,
      x = plotting_coords$mic_cut_x,
      xend = plotting_coords$mic_cut_xend,
      colour = "yellow", linewidth = 1.5
    ) +
    scale_y_discrete(limits = strain_order, labels = toupper(strain_order)) + # , guide = guide_axis(angle = 90)) +
    coord_equal()
}

# Bar plot scaled to match MIC heatmap
smg_plot <- function(smg_subset, strain_order) {
  smg_subset %>%
    filter(!strain %in% exclude) %>%
    arrange(factor(strain, levels = strain_order)) %>%
    ggplot(aes(y = strain, x = mean_48)) +
    geom_col(fill = "black") +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(family = "Helvetica", color = "black", size = 11),
      axis.text.x = element_text(family = "Helvetica", color = "black", size = 10),
      axis.line.y = element_line(color = "black", size = 0.2),
      axis.line.x = element_line(color = "black", size = 0.3),
      plot.title = element_text(family = "Helvetica", color = "black", size = 13)
    ) +
    xlab(NULL) +
    ylab(NULL) +
    # expand_limits(y = 1.0) +
    coord_fixed(ratio = 1.2) +
    scale_y_discrete(limits = strain_order) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "", "0.5", "", "1")
    )
}
