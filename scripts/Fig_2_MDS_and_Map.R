################################################################################
# Script to plot MDS plots of pair-wise distances
#
# Author:
#   Vladimir Bajić
#
# Figure(s):
#   - Figure 2
#   - Figure S2
#
################################################################################

# Packages ---------------------------------------------------------------------
library("tidyverse")
library("ape")
library("ggpubr")
library("rworldmap")
library("ggrepel")

# Set variables ----------------------------------------------------------------
path_fasta <- "out/fastas/1KGP_AFR_RSRS.fasta"
path_fullmetainfo <- "out/metadata_original_and_new_cluster_names.csv"
path_popadata <- "data/igsr_populations.tsv"
path_out_Fig2A <- "figures/Fig_2A_Map.pdf"
path_out_Fig2B <- "figures/Fig_2B_MDS.pdf"
path_out_Fig2 <- "figures/Fig_2_Map_and_MDS.pdf"
path_out_FigS2 <- "figures/Fig_S2_MDS.pdf"

# Load colors for plots --------------------------------------------------------
source("scripts/my_colors.R")
my_colors_gray <- my_colors
my_colors_gray[c(1:55)] <- "gray"

# Define shapes for plots ------------------------------------------------------
my_shape <-
  c(
  # African-American
  'ACB' = 2, 'ASW' = 6,
  # Eastern West Africa
  'MSL' = 0, 'GWD' = 5,
  # Western West Africa
  'ESN' = 3, 'YRI' = 4,
  # East Africa
  'LWK' = 1
  )

# Load ordered factors of cluster groups for plots -----------------------------
source("scripts/factor_lvls_ordered.R")

# Define ordering of clusterings for plots -------------------------------------
clustering_factor_lvls_ordered <-
  c(
    "SC", "SCL",
    "rhb_01", "rhb_02", "rhb_03",
    "tc_0.006", "tc_0.005", "tc_0.004", "tc_0.003"
  )

# Define function to plot MDS --------------------------------------------------
plot_mds <-
  function(mds, column) {
    p_mds <-
      mds %>%
      ggscatter(
        color = {{column}},
        title = column,
        x = "Dim1", y = "Dim2",
        label.rectangle = TRUE,
        palette = my_colors,
        shape = 20, size = 5, alpha = 1,
        ellipse = TRUE, ellipse.type = "convex"
      ) +
      theme(aspect.ratio = 1)
    return(p_mds)
  }

# Load function to plot hulls on ggplots ---------------------------------------
# Credit: https://stackoverflow.com/questions/18163153/convex-hulls-with-ggbiplot
source("scripts/func_StatBag.R")

# Load data --------------------------------------------------------------------
ape_fasta <- read.FASTA(path_fasta, type = "DNA")

# Load full metadata including secondary haplogroup groupings ------------------
metainfo <-
  read_csv(path_fullmetainfo) %>%
  select(!starts_with("original_")) %>%
  mutate_at(c("SCL", "SC"), as.factor) %>%
  filter(Sample != "RSRS") %>%
  mutate_if(is.numeric, as.factor) %>%
  as.data.frame()

# Update all factors with new unique levels ------------------------------------
metainfo[, c(-1, -2)] <-
  lapply(metainfo[, c(-1, -2)], factor, levels = factor_lvls_ordered)

# Calculate pairwise distances -------------------------------------------------
pwdist <-
  dist.dna(
    ape_fasta,
    model = "raw",
    variance = FALSE,
    gamma = FALSE,
    pairwise.deletion = FALSE,
    base.freq = NULL,
    as.matrix = TRUE
  )

# Calculate MDS ----------------------------------------------------------------
mds <-
  pwdist %>%
  cmdscale() %>%
  as_tibble() %>%
  mutate(Sample = rownames(pwdist)) %>%
  left_join(., metainfo, by = "Sample") %>%
  filter(Sample != "RSRS")

colnames(mds)[1:2] <- c("Dim1", "Dim2")

#-------------------------------------------------------------------------------

################################################################################
# Figure S2:  MDS plots for all groupings                                   ####
################################################################################

# MDS plots --------------------------------------------------------------------
mds_sc <- plot_mds(mds, "SC")
mds_scl <- plot_mds(mds, "SCL")
mds_rhb1 <- plot_mds(mds, "rhb_01")
mds_rhb2 <- plot_mds(mds, "rhb_02")
mds_rhb3 <- plot_mds(mds, "rhb_03")
mds_tc6 <- plot_mds(mds, "tc_0.006")
mds_tc5 <- plot_mds(mds, "tc_0.005")
mds_tc4 <- plot_mds(mds, "tc_0.004")
mds_tc3 <- plot_mds(mds, "tc_0.003")

# Combine all plots in one -----------------------------------------------------
p_all_mds_groupings <-
  ggarrange(
    ggpar(p = mds_sc, legend = "right"),
    ggpar(p = mds_scl, legend = "right"),
    ggpar(p = mds_rhb1, legend = "right"),
    ggpar(p = mds_rhb2, legend = "right"),
    ggpar(p = mds_rhb3, legend = "right"),
    ggpar(p = mds_tc6, legend = "right"),
    ggpar(p = mds_tc5, legend = "right"),
    ggpar(p = mds_tc4, legend = "right"),
    ggpar(p = mds_tc3, legend = "right"),
    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
    common.legend = FALSE, align = "hv"
  )

# Save plot --------------------------------------------------------------------
ggsave(path_out_FigS2, p_all_mds_groupings, width = 20, height = 20)


################################################################################
# Figure 2B:  MDS with SCL groupings                                        ####
################################################################################

p_mds <-
  mds %>%
  ggplot(., aes(Dim1, Dim2, colour = SCL, fill = SCL)) +
  theme_bw() +
  stat_bag(prop = 1, alpha = 0.2) +
  scale_fill_manual(values = my_colors_gray) +
  scale_colour_manual(values = my_colors_gray) +
  scale_shape_manual(values = my_shape) +
  theme(legend.position = "none") +
  annotate("text", x = -0.0015, y = -0.001, label = 'A, B, \nC, D, \nH, J, \nM, U', color = "gray", size = 8) +
  annotate("text", x = 0.002, y = -0.0005, label = 'L0', color = "gray", size = 8) +
  annotate("text", x = 0.0025, y = 0.0012, label = 'L1', color = "gray", size = 8) +
  annotate("text", x = -0.0005, y = 0.0015, label = 'L2', color = "gray", size = 8) +
  annotate("text", x = 0.000, y = -0.001, label = 'L3', color = "gray", size = 8) +
  annotate("text", x = -0.001, y = -0.0005, label = 'L4', color = "gray", size = 8) +
  annotate("text", x = 0.0005, y = 0.0002, label = 'L5', color = "gray", size = 8) +
  geom_point(aes(color = Population, shape = Population), size = 5, stroke = 2)

# Save plot --------------------------------------------------------------------
ggsave(path_out_Fig2B, p_mds, width = 8, height = 8)


################################################################################
# Figure 2A:  Population map                                                ####
################################################################################

# Get Map ----------------------------------------------------------------------
wmap <- getMap(resolution = "hight")

# Read in population metadata --------------------------------------------------
popadata <-
  path_popadata  %>%
  read_delim(., delim = "\t", show_col_types = FALSE) %>%
  filter(`Superpopulation code` == "AFR") %>%
  select(
    "Population code",
    "Superpopulation code",
    "Population latitude",
    "Population longitude"
  )

# Plot Map ---------------------------------------------------------------------
p_map <-
  ggplot() +
  geom_polygon(
    data = wmap,
    aes(x = long, y = lat, group = group),
    fill = "lightgray", color = "white"
  ) +
  coord_cartesian(xlim = c(-100, 40), ylim = c(-5, 40)) +
  geom_point(
    data = popadata,
    aes(
      x = `Population longitude`,
      y = `Population latitude`,
      color = `Population code`,
      shape = `Population code`
    ),
    size = 5, stroke = 2
  ) +
  geom_label_repel(
    data = popadata,
    aes(
      x = `Population longitude`,
      y = `Population latitude`,
      color = `Population code`,
      label = `Population code`
    ),
    size = 8
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = my_colors_gray) +
  scale_colour_manual(values = my_colors_gray) +
  scale_shape_manual(values = my_shape)

# Save plot --------------------------------------------------------------------
ggsave(path_out_Fig2A, p_map, width = 10, height = 5)


################################################################################
# Figure 2:  Panel A: Map | Panel B: MDS plot                               ####
################################################################################

# Combine Map and MDS in one plot ----------------------------------------------
p_all <-
  ggarrange(
    p_map, p_mds,
    labels = c("A", "B"),
    common.legend = FALSE, align = "v",
    ncol = 1, legend = "none", heights = c(1, 2)
  )

# Save combined plot -----------------------------------------------------------
ggsave(path_out_Fig2, p_all, width = 10, height = 14.8)