################################################################################
# Script to make CA plots
#
# Author:
#   Vladimir Bajić
#
# Figure(s):
#   - Figure 6
#
################################################################################

# Packages ---------------------------------------------------------------------
library("tidyverse")
library("FactoMineR")
library("factoextra")
library("ggpubr")

# Set variables ----------------------------------------------------------------
path_all_metadata <- "out/metadata_original_and_new_cluster_names.csv"
path_out <- "figures/Fig_6_CA.pdf"

# Define function to plot CA plots ---------------------------------------------
plot_ca <-
  function(meta_data, column, pop_colors) {
    ca <- meta_data %>%
    select(Population, {{column}}) %>%
    droplevels() %>%
    table() %>%
    prop.table(margin = 1) %>%
    CA(graph = FALSE)  %>%
    fviz_ca_biplot(
      col.col = "grey", 
      col.row = as.list(unname(pop_colors)),
      labelsize = 5,
      pointsize = 5,
      arrow = c(FALSE, TRUE),
      repel = TRUE, title = column
      ) +
    theme(aspect.ratio = 1, legend.position = "none")
    return(ca)
  }

# Load colors for plots --------------------------------------------------------
source("scripts/my_colors.R")
pop_colors <- my_colors[56:62] # get only populations

# Load ordered factors of cluster groups for plots -----------------------------
source("scripts/factor_lvls_ordered.R")

# Set ggrepel options ----------------------------------------------------------
# increase accepted overlap for plotting
options(ggrepel.max.overlaps = Inf)

# Load all metadata  -----------------------------------------------------------
# And update all factors with new levels
all_metadata <-
  read_csv(path_all_metadata) %>%
  filter(Sample != "RSRS")  %>%
  mutate_if(is.numeric, as.factor) %>%
  mutate(across(all_of(c(starts_with("rhb_"), starts_with("tc_"), starts_with("SC"))), ~factor(.x, levels = factor_lvls_ordered))) %>%
  mutate(Population = factor(Population, levels = c("ASW", "ACB", "GWD", "MSL", "YRI", "ESN", "LWK"))) %>%
  as.data.frame()

# CA plots ---------------------------------------------------------------------
ca_sc <- plot_ca(all_metadata, "SC", pop_colors)
ca_scl <- plot_ca(all_metadata, "SCL", pop_colors)
ca_tc_3 <- plot_ca(all_metadata, "tc_0.003", pop_colors)
ca_tc_4 <- plot_ca(all_metadata, "tc_0.004", pop_colors)
ca_tc_5 <- plot_ca(all_metadata, "tc_0.005", pop_colors)
ca_tc_6 <- plot_ca(all_metadata, "tc_0.006", pop_colors)
ca_rhb_02 <- plot_ca(all_metadata, "rhb_02", pop_colors)
ca_rhb_03 <- plot_ca(all_metadata, "rhb_03", pop_colors)

# Combine all plots in one -----------------------------------------------------
ca_all <-
  ggarrange(
    ca_sc, ca_scl,
    ca_rhb_02, ca_rhb_03,
    ca_tc_6, ca_tc_5,
    ca_tc_4, ca_tc_3,
    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
    common.legend = FALSE, legend = "none", ncol = 2, nrow = 4
  )

# Save plot --------------------------------------------------------------------
ggsave(path_out, ca_all, width = 10, height = 20)