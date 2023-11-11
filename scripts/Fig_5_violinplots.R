################################################################################
# Script to plot violin plots of pair-wise distances
#
# Author:
#   Vladimir Bajić
#
# Figure(s):
#   - Figure 5
#
################################################################################

# Packages ---------------------------------------------------------------------
library("tidyverse")
library("ape")
library("ggpubr")

# Set variables ----------------------------------------------------------------
path_fasta <- "out/fastas/1KGP_AFR_RSRS.fasta"
path_fullmetainfo <- "out/metadata_original_and_new_cluster_names.csv"
path_out <- "figures/Fig_5_Violin.pdf"

# Load colors for plots --------------------------------------------------------
source("scripts/my_colors.R")

# Load ordered factors of cluster groups for plots -----------------------------
source("scripts/factor_lvls_ordered.R")

# Define ordering of clusterings for plots -------------------------------------
clustering_factor_lvls_ordered <-
  c(
    "SC", "SCL",
    "rhb_01", "rhb_02", "rhb_03",
    "tc_0.006", "tc_0.005", "tc_0.004", "tc_0.003"
  )

# Load data --------------------------------------------------------------------
ape_fasta <- read.FASTA(path_fasta, type = "DNA")

# Load full metadata including secondary haplogroup groupings ------------------
metainfo <-
  read_csv(path_fullmetainfo) %>%
  select(
    Sample,
    rhb_01, rhb_02, rhb_03,
    tc_0.006, tc_0.005, tc_0.004, tc_0.003,
    SC, SCL
  ) %>%
  mutate_at(c("SCL", "SC"), as.factor) %>%
  filter(Sample != "RSRS") %>%
  mutate_if(is.numeric, as.factor) %>%
  as.data.frame()

# Update all factors with new unique levels ------------------------------------
metainfo[, -1] <- lapply(metainfo[, -1], factor, levels = factor_lvls_ordered)

# Calculate pairwise distances -------------------------------------------------
res_table <-
  dist.dna(
    ape_fasta,
    model = "raw",
    variance = FALSE,
    gamma = FALSE,
    pairwise.deletion = FALSE,
    base.freq = NULL,
    as.matrix = TRUE
  ) %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  gather("Var2", "value", -Var1) %>%
  filter(Var1 != "RSRS" & Var2 != "RSRS")

# Add to the results table of pairwise distances the info
# about which groupings each part of the pair is in
res_table_groupings <-
  left_join(res_table, metainfo, by = c("Var1" = "Sample")) %>%
  left_join(., metainfo, by = c("Var2" = "Sample")) %>%
  pivot_longer(cols = ends_with(".x"),
               names_to = "clustering.x", values_to = "cluster.x") %>%
  pivot_longer(cols = ends_with(".y"), 
               names_to = "clustering.y", values_to = "cluster.y") %>%
  mutate(clustering.x = str_remove(clustering.x, ".x"),
         clustering.y = str_remove(clustering.y, ".y"))

# Select only those individuals that belong to the same cluster ----------------
vptmp <-
  res_table_groupings %>%
  filter(clustering.x == clustering.y) %>%
  unite("cluster_names.x", c(clustering.x, cluster.x),
        remove = FALSE, sep = "_") %>%
  unite("cluster_names.y", c(clustering.y, cluster.y),
        remove = FALSE, sep = "_") %>%
  filter(cluster_names.x == cluster_names.y) %>%
  mutate(clustering.x = as.factor(clustering.x)) %>%
  # Update levels (order of clustering for plot)
  mutate(clustering.x = factor(clustering.x, levels = clustering_factor_lvls_ordered),
         cluster.x = factor(cluster.x, levels = factor_lvls_ordered))

# Make violin plot of pair wise differences ------------------------------------
p_violin_facet <-
  vptmp %>%
  ggplot(aes(x = reorder(cluster.x, value), fill = cluster.x, y = value)) +
  theme_bw() +
  coord_flip() +
  geom_violin() +
  scale_fill_manual(values = my_colors) +
  geom_point(stat = "summary", fun = "mean") +
  facet_grid(clustering.x ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ylab("pairwise distance") +
  xlab(element_blank()) +
  geom_hline(yintercept = c(0.003, 0.004, 0.005, 0.006), linetype = "dashed" , alpha = 0.5) +
  theme(strip.placement = "outside", legend.position = "none", strip.text.y.left = element_text(angle = 0))

# Save plot --------------------------------------------------------------------
ggsave(path_out, p_violin_facet, width = 5, height = 15)