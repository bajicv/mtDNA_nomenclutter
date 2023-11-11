################################################################################
# Script to plot barplots of haplogroup frequencies across diffrent groupings
#
# Author:
#   Vladimir Bajić
#
# Figure(s):
#   - Figure 4
#
################################################################################

# Packages ---------------------------------------------------------------------
library("tidyverse")
library("cowplot")
library("ggpubr")

# Set variables ----------------------------------------------------------------
path_all_metadata <- "out/metadata_original_and_new_cluster_names.csv"
path_out <- "figures/Fig_4_Barplot.pdf"

# Define function to plot barplots ---------------------------------------------
plot_barplot <-
  function(meta_data, column) {
    bp_0 <-
      meta_data  %>%
      group_by(Population)  %>%
      ggplot(aes(Population, fill = as.factor(meta_data[, column]))) +
      geom_bar(position = "fill") +
      theme_minimal() +
      scale_fill_manual(values = my_colors, na.value = "white") +
      coord_flip() +
      labs(title = as.character(column), fill = as.character(column), y = "Frequency", x = NULL)

  bp_1 <- bp_0 + theme(legend.position = "none")
  bp_legend <- get_legend(bp_0 + theme(legend.box.margin = margin(0, 0, 0, 12)))
  bp <- plot_grid(bp_1, bp_legend, rel_widths = c(4, 1))
  return(bp)
}

# Load colors for plots --------------------------------------------------------
source("scripts/my_colors.R")

# Load ordered factors of cluster groups for plots -----------------------------
source("scripts/factor_lvls_ordered.R")

# Load all metadata  -----------------------------------------------------------
# And update all factors with new levels
all_metadata <-
  read_csv(path_all_metadata) %>%
  filter(Sample != "RSRS")  %>%
  mutate_if(is.numeric, as.factor) %>%
  mutate(across(all_of(c(starts_with("rhb_"), starts_with("tc_"), starts_with("SC"))), ~factor(.x, levels = factor_lvls_ordered))) %>%
  mutate(Population = factor(Population, levels = c("ASW", "ACB", "GWD", "MSL", "YRI", "ESN", "LWK", "mt-MRCA"))) %>%
  as.data.frame()

# Bar plots --------------------------------------------------------------------
bp_sc <- plot_barplot(all_metadata, "SC")
bp_scl <- plot_barplot(all_metadata, "SCL")
bp_rhb_03 <- plot_barplot(all_metadata, "rhb_03")
bp_rhb_02 <- plot_barplot(all_metadata, "rhb_02")
bp_rhb_01 <- plot_barplot(all_metadata, "rhb_01")
bp_tc_6 <- plot_barplot(all_metadata, "tc_0.006")
bp_tc_5 <- plot_barplot(all_metadata, "tc_0.005")
bp_tc_4 <- plot_barplot(all_metadata, "tc_0.004")
bp_tc_3 <- plot_barplot(all_metadata, "tc_0.003")

# Combine all plots in one -----------------------------------------------------
bp_all <-
  ggarrange(
    bp_sc, bp_scl, bp_rhb_01,
    bp_rhb_02, bp_rhb_03, bp_tc_6,
    bp_tc_5, bp_tc_4, bp_tc_3,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    common.legend = FALSE,
    legend = "none"
  )

# Save plot --------------------------------------------------------------------
ggsave(path_out, bp_all, width = 20, height = 20)