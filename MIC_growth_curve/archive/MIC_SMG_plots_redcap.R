## ---------------------------
## Purpose: Calculate, plot and upload MIC and SMG data to REDCap
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## OD readings are read from one sheet and pivoted to "long" format,
## metadata in long format is read from another sheet and concatenated
## ---------------------------
options(scipen = 999)
## ---------------------------
# load packages
library(jsonlite)
## ---------------------------

# Functions
source("MIC_growth_curve/MIC_calc_functions.R")
source("MIC_growth_curve/MIC_heatmap.R")

# VARIABLES - READ THESE CAREFULLY!!!!!!!!!!!!!!!!!!!!!!
# For uploading to REDcap
api_url <- "https://redcap.ahc.umn.edu/redcap/api/"

# Enter drug and spreadsheet data
input_drug <- "flc"
replicates <- 3

mic_spreadsheet <- ""
smg_spreadsheet <- ""

od_tab <- 1 # Excel tab of OD values
metadata_tab <- 2 # Excel tab of metadata
metadata_range <- "A1:F13"

# Type either "strain" or "concentration" for col names, depending on plate layout
column_names <- "concentration"

# Number of columns filled, default is 12
cols_used <- 12

# Drug concentrations (this assumes the usual "screening" set-up, CHANGE AS NEEDED)
concentration <- case_when(
  toupper(input_drug) == "FLC" ~ c(0, 0.5, 1, 2, 4, 8, 16, 32),
  toupper(input_drug) %in% c("MCF", "AMB") ~ c(0, 0.016, 0.032, 0.064, 0.125, 0.256, 0.5, 1)
)

max_concentration <- 32

# Add control IDs
control_strains <- c("AMS5123", "AMS5122")

# FOR SPECIFIC ORDER: put IDs in "strains" vector. Otherwise leave commented out.
# strain_list <- c("AMS5123", "MEC257", "MEC208", "MEC205", "MEC194", "MEC190", "MEC188", "MEC176", "MEC175")

# TO SKIP STRAINS IN PLOT:
# Need strain_list above with wanted strains plus IDs to skip in "exclude" vector. Otherwise leave commented out.
# exclude <- c("MEC174", "MEC185", "MEC196")

mic_date <- str_extract(mic_spreadsheet, "\\d+-\\d+-\\d+")

################################################################################
# Load metadata
meta.frame <- read_excel(mic_spreadsheet,
  sheet = metadata_tab,
  range = metadata_range
)

meta.frame$drug <- c(toupper(input_drug), rep(NA, times = cols_used - length(input_drug)))

meta.frame$concentration <- c(concentration, rep(NA, times = cols_used - length(concentration)))

# Set vars for use below
meta_names <- case_when(
  column_names == "strain" ~ "strain",
  column_names == "concentration" ~ "concentration"
)
meta_col <- case_when(
  column_names != "strain" ~ "strain",
  column_names != "concentration" ~ "concentration"
)

cutoff <- case_when(
  meta.frame$drug[1] %in% c("FLC", "MCF") ~ 0.5,
  meta.frame$drug[1] == "AMB" ~ 0.1
)

# x_axis_angle <- case_when(meta.frame$drug[1]== "FLC" ~ 0,
#                          meta.frame$drug[1] == "MCF" ~ 90,
#                          meta.frame$drug[1] == "AMB" ~ 90)
x_axis_angle <- 0

################################################################################
# MIC
# Loop over each plate
for (i in 1:replicates) {
  d.frame <- read_excel(mic_spreadsheet,
    sheet = od_tab,
    range = meta.frame[[tolower(input_drug)]][i]
  )

  # Colnames set from input variable
  colnames(d.frame) <- as.character(pull(meta.frame, column_names))

  # Add (strain or concentration) not used as column names as a new column
  if (!all(meta.frame$strain[!is.na(meta.frame$strain)] %in% colnames(d.frame))) {
    d.frame$strain <- meta.frame$strain[!is.na(meta.frame$strain)]
  } else {
    d.frame$concentration <- as.character(meta.frame$concentration[!is.na(meta.frame$concentration)])
  }

  # Add remaining metadata
  d.frame <- d.frame %>%
    add_column(plate = i) %>%
    add_column(drug = meta.frame$drug[!is.na(meta.frame$drug)]) %>%
    add_column(media = meta.frame$media[!is.na(meta.frame$media)]) %>%
    add_column(temp = meta.frame$temp[!is.na(meta.frame$temp)])

  assign(paste0("rep", i), d.frame)
}

# Pull together and tidy data
drug_data <- bind_rows(mget(ls(pattern = "^rep\\d+$")))
drug <- pivot_longer(drug_data,
  names_to = meta_names,
  values_to = "OD600",
  cols = -c(all_of(meta_col), plate, drug, media, temp)
)

# BLANKS: this checks for "blank" labeled wells in metadata
# if none found, sets it to the highest concentration well of the control strain
drug <- drug %>%
  mutate(concentration = case_when((!("blank" %in% c(strain, concentration)) & concentration == as.character(max_concentration) & strain %in% control_strains) ~ "blank",
    .default = concentration
  ))
drug$concentration <- factor(tolower(drug$concentration), levels = mixedsort(unique(drug$concentration)))

# Create strain list for plotting (alphanumeric + control at bottom) unless it's been entered above
if (exists("strain_list")) {
  strains <- strain_list
} else {
  strains <- c(unique(drug$strain[drug$strain %in% control_strains]), rev(unique(drug$strain[!(drug$strain %in% control_strains)])))
}

mic_boxplot(drug)
################################################################################
# Read in SMG ODs and reuse metadata from MIC spreadsheet
smg_range <- read_excel(smg_spreadsheet, sheet = metadata_tab)
for (i in 1:replicates) {
  d.frame <- read_excel(smg_spreadsheet,
    sheet = od_tab,
    range = smg_range[[tolower(input_drug)]][i]
  )

  # Colnames set from input variable
  colnames(d.frame) <- as.character(pull(meta.frame, column_names))

  # Add (strain or concentration) not used as column names as a new column
  if (!all(colnames(d.frame) == meta.frame$strain)) {
    d.frame$strain <- meta.frame$strain[!is.na(meta.frame$strain)]
  } else {
    d.frame$concentration <- as.character(meta.frame$concentration[!is.na(meta.frame$concentration)])
  }

  # Add remaining metadata
  d.frame <- d.frame %>%
    add_column(plate = i) %>%
    add_column(drug = meta.frame$drug[!is.na(meta.frame$drug)]) %>%
    add_column(media = meta.frame$media[!is.na(meta.frame$media)]) %>%
    add_column(temp = meta.frame$temp[!is.na(meta.frame$temp)])

  assign(paste("rep", i, sep = ""), d.frame)
}

# Pull together and tidy data
drug48_data <- bind_rows(mget(ls(pattern = "^rep\\d+$")))
drug48 <- pivot_longer(drug48_data,
  names_to = meta_names,
  values_to = "OD600",
  cols = -c(all_of(meta_col), plate, drug, media, temp)
)

# This sets the blank when it's one well of the control strain
drug48 <- drug48 %>%
  mutate(concentration = case_when((!("blank" %in% c(strain, concentration)) & concentration == as.character(max_concentration) & strain %in% control_strains) ~ "blank",
    .default = concentration
  ))
drug48$concentration <- factor(tolower(drug48$concentration), levels = mixedsort(unique(drug48$concentration)))

mic_boxplot(drug48)
################################################################################
# Calculate relative ODs and MIC values
drug_od <- calculate_od(drug)

drug_mic <- mic(drug_od, cutoff)
drug_smg_input <- smg_input(drug_od, cutoff)

drug_smg_od <- calculate_od(drug48)
drug_smg <- smg_subset(drug_smg_od, drug_smg_input)

################################################################################
# Finally plot

# Get the points for yellow MIC line
drug_plotting_coords <- plotting_coords(drug_od, drug_mic, cutoff, strains)

# Make the "heatmap"
drug_mic_plot <- mic_plot(drug_od, strains, drug_plotting_coords) +
  # xlab(case_when(drug$drug[1] %in% c("AMB", "MCF") ~ paste("\n", drug$drug[1]," MIC"),
  #    .default = paste(drug$drug[1], " MIC"))) +
  xlab(paste(drug$drug[1], " MIC")) +
  theme(
    axis.text.x = element_text(angle = x_axis_angle),
    plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"),
    axis.title.y = element_text(
      hjust = 0.5,
      margin = margin(0, 10, 0, 0)
    )
  )

# Barplot
drug_smg_plot <- smg_plot(left_join(drug_mic, drug_smg), strains) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(),
    plot.margin = unit(c(0, 0.5, 0.5, 0), "cm")
  ) + # ,
  # axis.title.x = element_text(lineheight = 1.1)) +
  # xlab(case_when(drug$drug[1] %in% c("AMB", "MCF")~ paste("\n\n","SMG"),
  # .default = "SMG"))
  xlab("SMG")

# Pull together with patchwork
drug_full_plot <- (drug_mic_plot + drug_smg_plot) + plot_layout(guides = "collect")

# Save as desired
ggsave(paste0("images/2023_MICs/", mic_date, "_MEC_", drug$drug[1], "_MIC.png"),
  drug_full_plot,
  width = 7.3,
  height = 8,
  units = "in",
  device = png,
  bg = "white",
  dpi = 300
)
################################################################################
# REDCap imports
# Numeric values to text

#### FILL IN QC RESULTS HERE #################
redcap_qc <- "" # 1 for yes, 0 for no

drug_mic$concentration <- as.character(drug_mic$concentration)

# Append > when needed
drug_mic <- drug_mic %>%
  mutate(concentration = case_when((mean_norm_OD > cutoff & concentration == max_concentration) ~ paste0(">", max_concentration),
    .default = concentration
  ))

redcap_drug_value <- case_when(
  drug$drug[1] == "FLC" ~ 0,
  drug$drug[1] == "MCF" ~ 1,
  drug$drug[1] == "AMB" ~ 2
)

redcap_drug_man <- case_when(
  drug$drug[1] == "FLC" ~ "Alfa Aesar",
  drug$drug[1] == "MCF" ~ "MedChemExpress",
  drug$drug[1] == "AMB" ~ "Chem-Impex Intl Inc"
)

redcap_drug_lot <- case_when(
  drug$drug[1] == "FLC" ~ "P04F029",
  drug$drug[1] == "MCF" ~ "22367",
  drug$drug[1] == "AMB" ~ "001448-180229"
)

redcap_drug_solvent <- case_when(
  drug$drug[1] == "FLC" ~ "EtOH",
  drug$drug[1] == "MCF" ~ "DMSO",
  drug$drug[1] == "AMB" ~ "DMSO"
)

# For each strain, create new record and send form
for (i in c(1:length(drug_mic$strain))) {
  primary_id <- drug_mic$strain[i]

  record <- c(
    primary_id = primary_id,
    redcap_repeat_instrument = "mic_results",
    redcap_repeat_instance = "new",
    drug = redcap_drug_value,
    mic_date = mic_date,
    mic50 = as.character(drug_mic$concentration[drug_mic$strain == primary_id]),
    mic_file_name = basename(mic_spreadsheet),
    smg = as.character(drug_smg$mean_48[drug_smg$strain == primary_id]),
    smg_file_name = basename(smg_spreadsheet),
    mic_media = 1,
    mic_temp = 1,
    qc_ok = redcap_qc,
    drug_source = redcap_drug_man,
    lot_number = redcap_drug_lot,
    solvent = redcap_drug_solvent,
    mic_results_complete = 2
  )

  result_data <- toJSON(list(as.list(record)), auto_unbox = TRUE)

  formData <- list(
    "token" = Sys.getenv("redcap_api_key"),
    content = "record",
    action = "import",
    format = "json",
    type = "flat",
    overwriteBehavior = "normal",
    forceAutoNumber = "false",
    data = result_data,
    returnContent = "count",
    returnFormat = "json"
  )
  response <- httr::POST(api_url, body = formData, encode = "form")
  result <- httr::content(response)
  print(result)
}
