################################################################################
# Script to plot tree with comparison between
# nomenclature-based (NBG) and algorithm-based (ABG) groupings.
#
# Author:
#   Vladimir Bajić
#
# Figure(s):
#   - Figure 3
#   - Figure S1
#
################################################################################

# Packages ---------------------------------------------------------------------
library(tidyverse)
library(ggtree)
library(ape)

# Set variables ----------------------------------------------------------------
path_tree <- "out/trees/1KGP_AFR_RSRS_iqtree.contree"
path_fullmetainfo <- "out/metadata_original_and_new_cluster_names.csv"
path_tree_summary_metadata <- "out/metadata_tree_summary.csv"
path_out_FigS1 <- "figures/Fig_S1_Tree.pdf"
path_out_Fig3 <- "figures/Fig_3_Tree.pdf"

# Load colors for plots --------------------------------------------------------
source("scripts/my_colors.R")

# Load ordered factors of cluster groups for plots -----------------------------
source("scripts/factor_lvls_ordered.R")

# # Load metadata --------------------------------------------------------------
metadata <-
  read_csv(path_fullmetainfo) %>%
  select(!starts_with("original_")) %>%
  mutate_at(c("SCL", "SC"), as.factor) %>%
  filter(Sample != "RSRS") %>%
  mutate_if(is.numeric, as.factor) %>%
  as.data.frame()

# Load tree summary metadata ---------------------------------------------------
tree_summary_metadata <-
  read_delim(path_tree_summary_metadata) %>%
  mutate_at(c("SCL", "SC"), as.factor)

# Load iqtree and add metadata -------------------------------------------------
iqtree <-
  read.tree(path_tree) %>%
  left_join(., metadata, by = c("label" = "Sample"))

# Set RSRS Haplogorups as mt-MRCA
iqtree@extraInfo[1, "Haplogroup"] <- "mt-MRCA"

################################################################################
# Figure S1:  Full AFR Tree                                                 ####
################################################################################

# Define which nodes belong to which SCL group (needed for groupOTU)
groupInfo <- split(iqtree@extraInfo$node, iqtree@extraInfo$SCL)

# First group OTU based on SCL groups
# This gorup info is then used for coloring "color=group"
p_FigS1 <-
  iqtree %>%
  groupOTU(., groupInfo) %>%
  left_join(., metadata, by = c("label" = "Sample")) %>%
  ggtree(layout = "fan", aes(color=group)) +
  geom_tiplab(aes(label = Haplogroup), align = TRUE, size = 2.5) +
  scale_color_manual(values = my_colors, na.value = "black") +
  theme(legend.position = "none")

# Save plot --------------------------------------------------------------------
ggsave(path_out_FigS1, p_FigS1, width = 20, height = 20)

#-------------------------------------------------------------------------------


################################################################################
# Figure 3:   Reduced phylogeny                                             ####
################################################################################

# Make reduced iqtree ----------------------------------------------------------

# Make list of all individuals to keep in the tree
tokeep <- c("RSRS", tree_summary_metadata$Sample)

# Subset tree
iqtree_reduced <-
  iqtree@phylo %>%
  keep.tip(., tokeep)  %>%
  left_join(., tree_summary_metadata, by = c("label" = "Sample"))

# Set RSRS Haplogorups and short_haplogroup as mt-MRCA
iqtree_reduced@extraInfo[1, "Haplogroup"] <- "mt-MRCA"
iqtree_reduced@extraInfo[1, "short_haplogroup"] <- "mt-MRCA"

# Save reduced tree
gg_red <-
  iqtree_reduced  %>%
  ggtree(layout = "rectangular") +
  geom_tiplab(aes(label = short_haplogroup), align = TRUE)

# Check if some nodes should be fliped
# gg_red + geom_text(aes(label = node))

# Flip branches so that non-L haplogroups are at the top
gg_red_flip <- gg_red  %>%
  flip(52, 48)  %>%
  flip(55, 53)  %>%
  flip(23, 24)  %>%
  flip(64, 62)

#-------------------------------------------------------------------------------

# Extract columns that are relevant for plot
# And update all factors with factor_lvls_ordered
toplot <-
  tree_summary_metadata  %>%
  select(
    Sample, SC, SCL,
    `rhb_01`, `rhb_02`, `rhb_03`,
    `tc_0.006`, `tc_0.005`, `tc_0.004`, `tc_0.003`
  )  %>%
  mutate_if(is.numeric, as.factor)  %>%
  column_to_rownames("Sample")  %>%
  mutate(across(everything(), ~factor(.x, levels = factor_lvls_ordered)))

# Calculate number of individuals in the dataset
# with the same combination of clusterings
n_ind_per_hap <-
  tree_summary_metadata  %>%
  select(Sample, starts_with("n_"))  %>%
  pivot_longer(cols = starts_with("n_")) %>%
  group_by(Sample)  %>%
  summarise(value = sum(value))  %>%
  ungroup()

# Plot reduced tree and heatmap
# with number of individuals on the right side of the tree
p_Fig3_tmp <-
  gheatmap(
    gg_red_flip,
    toplot,
    offset = 0.001,
    width = 1,
    colnames = TRUE,
    colnames_angle = -90
  ) +
  scale_fill_manual(values = my_colors, na.value = "white") +
  theme(legend.position = "none") +
  geom_facet(
    panel = "N = 660", data = n_ind_per_hap, geom = geom_col, aes(x = n),
    fill = "gray30", color = "gray30", orientation = "y", width = 0.8
  ) +
  geom_facet(
    panel = "N = 660", data = n_ind_per_hap, geom = geom_text,
    aes(x = 0 - 20, label = as.character(value))
  )

# Define how much space each of elements gets
p_Fig3 <- facet_widths(p_Fig3_tmp, c(Tree = 4))

# Save plot --------------------------------------------------------------------
ggsave(path_out_Fig3, p_Fig3, width = 15, height = 20)