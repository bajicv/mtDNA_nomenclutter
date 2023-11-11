################################################################################
# Script to create metadata for the mtDNA from AFR superpop from 1KGP
#
# This script:
#   - removes known related individuals
#   - adds HaploGrep3 results
#   - adds rhierBAPS results
#   - adds TreeCluster results
#
# Author:
#   Vladimir Bajić
#
################################################################################

# Packages ---------------------------------------------------------------------
library(tidyverse)

# Set variables ----------------------------------------------------------------
path_out_table <- "out/metadata.txt"
path_out_tree_summary_metadata <- "out/metadata_tree_summary.csv"
path_out_all_metadata <- "out/metadata_original_and_new_cluster_names.csv"
path_popmetadata <- "data/igsr_samples.tsv"
path_relatives <- "data/20140625_related_individuals.txt"
path_treecluster_out_dir <- "out/treecluster"
path_rhierBAPS_out <- "out/rhierBAPS/rhb_1KGP_AFR_level10_keep_singletons_T.csv"
path_haplogrep3_out <- "out/haplogrep3/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_noGaps-haplogrep3_phylotree-fu-rcrs-1.2.txt"

# Define functions  ------------------------------------------------------------

## Load function to find longest shared name of haplogorups
source("scripts/func_lcPrefix.R")

## func to match clusters
source("scripts/func_match_clusters.R")

## Define function to rename -1 TreeClusters
tc_singleton_rename <-
  function(df) {
    df[df$ClusterNumber == -1, 2] <- seq(
      max(df$ClusterNumber) + 1,
      length = nrow(df[df$ClusterNumber == -1, ])
    )
    return(df)
  }

## func to left join columns based on key-value data
# used to create full metadata with renamed clusters
ljc <-
  function(full_data, key_val_data, column_to_join) {
    all_data <- left_join(
      full_data,
      unique(
        key_val_data[, c(column_to_join, paste0("original_", column_to_join))]
      )
    )
    return(all_data)
  }

# Load metadata ----------------------------------------------------------------
# downloaded from: https://www.internationalgenome.org/data-portal/sample
igsr_samples <-
  read_tsv(path_popmetadata) %>%
  rename(Sample =  "Sample name", Population = "Population code") %>%
  select(Sample, Population)

# Load known relatives ---------------------------------------------------------
# downloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
relatives <- read_delim(path_relatives)

# Load HaploGrep3 output -------------------------------------------------------
# from running HaploGrep3 on
# 1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_noGaps.fa
haplogrep <-
  read_delim(path_haplogrep3_out) %>%
  rename(Sample = "SampleID") %>%
  select(Sample, Haplogroup)

# Load TreeCluster output ------------------------------------------------------
# from running TreeCluster on 1KGP_fasta_RSRS_african_iqtree_tree.contree
# and rename -1 TreeClusters
file_paths <-
  list.files(
    path = path_treecluster_out_dir,
    pattern = "*.txt",
    full.names = TRUE
  )

tc_6 <-
  read_delim(file_paths[grep("t_0.006.txt", file_paths)]) %>%
  tc_singleton_rename()

tc_5 <-
  read_delim(file_paths[grep("t_0.005.txt", file_paths)]) %>%
  tc_singleton_rename()

tc_4 <-
  read_delim(file_paths[grep("t_0.004.txt", file_paths)]) %>%
  tc_singleton_rename()

tc_3 <-
  read_delim(file_paths[grep("t_0.003.txt", file_paths)]) %>%
  tc_singleton_rename()

# Left join all of them inot one table
tc <-
  reduce(
    list(tc_6, tc_5, tc_4, tc_3),
    left_join, by = "SequenceName"
  ) %>%
  rename_with(
    ~ c("Sample", "tc_0.006",  "tc_0.005",  "tc_0.004",  "tc_0.003"),
    everything()
  )

# Load rhierBAPS output --------------------------------------------------------
rhb <-
  read_delim(path_rhierBAPS_out) %>%
  select(Isolate, `level 1`, `level 2`, `level 3`) %>%
  rename(Sample = "Isolate", rhb_01 = "level 1", rhb_02 = "level 2", rhb_03 = "level 3")

# Secondary nomenclature-based haplogorup groupings (NBG) ####
# create SC groups by extracting first letter of each haplogroup
# create SCL groups by extracting first 2 letters for L haplogroups and only first letter for other haplogroups 
nbg <-
  haplogrep %>%
  mutate("SC" = substr(Haplogroup, 1, 1)) %>%
  mutate(
    "SCL" = case_when(startsWith(Haplogroup, "L") ~ substr(Haplogroup, 1, 2),
    !startsWith(Haplogroup, "L") ~ substr(Haplogroup, 1, 1))
  ) %>%
  select(Sample, SC, SCL)

# Mergeing and filtering for relatives
metadata <-
  left_join(tc, rhb) %>%
  left_join(., nbg) %>%
  left_join(., haplogrep) %>%
  left_join(., igsr_samples) %>%
  filter(!(Sample %in% relatives$Sample)) %>%
  select(
    Sample, Population, Haplogroup,
    SC, SCL, rhb_01, rhb_02, rhb_03,
    tc_0.006, tc_0.005, tc_0.004, tc_0.003
  )

metadata[metadata$Sample == "RSRS", "Population"] <- "mt-MRCA"
metadata[metadata$Sample == "RSRS", "Haplogroup"] <- "mt-MRCA"
metadata[metadata$Sample == "RSRS", "SC"] <- "mt-MRCA"
metadata[metadata$Sample == "RSRS", "SCL"] <- "mt-MRCA"

# Save metadata table ----------------------------------------------------------
# write_delim(metadata, path_out_table)

################################################################################

# Make summary for all results -------------------------------------------------
summary <-
  metadata %>%
  group_by(across(c(-Sample, -Population, -Haplogroup))) %>%
  summarise(
    n = n(),
    Samples = list(Sample),
    POPs = list(Population),
    Haplogroups = list(Haplogroup),
    short_haplogroup = lcPrefix(unlist(Haplogroup)),
    n_mtMRCA = sum(unlist(POPs) == "mt-MRCA"),
    n_ASW = sum(unlist(POPs) == "ASW"),
    n_ESN = sum(unlist(POPs) == "ESN"),
    n_ACB = sum(unlist(POPs) == "ACB"),
    n_LWK = sum(unlist(POPs) == "LWK"),
    n_MSL = sum(unlist(POPs) == "MSL"),
    n_GWD = sum(unlist(POPs) == "GWD"),
    n_YRI = sum(unlist(POPs) == "YRI")
  ) %>%
  ungroup()

names(metadata)

# Manuall haplogroup name adjustments ------------------------------------------
# Change ambiguous haplogroups after inspecting phylogenies at Haplogrep3

## D1*4 -> simplify the name
summary$Haplogroups[4]
summary$short_haplogroup[4]
summary$short_haplogroup[4] <- "D1"

## L0b* -> simplify the name
summary$Haplogroups[8]
summary$short_haplogroup[8]
summary$short_haplogroup[8] <- "L0b"

## 1xL1c1b 4xL1c1a 6xL1c1d -> L1c1a'b'd
summary$Haplogroups[14]
summary$short_haplogroup[14]
summary$short_haplogroup[14] <- "L1c1a'b'd"

## L2a5*1 -> simplify the name
summary$Haplogroups[20]
summary$short_haplogroup[20]
summary$short_haplogroup[20] <- "L2a5"

## 3xL3k vs 1xL3i -> take L3k
summary$Haplogroups[24]
summary$short_haplogroup[24]
summary$short_haplogroup[24] <- "L3k"

## 2xL3x vs 98xL3e -> take L3e
summary$Haplogroups[28]
summary$short_haplogroup[28]
summary$short_haplogroup[28] <- "L3e"


# Data preparation for plotting ------------------------------------------------

# Sample only one individual from all groups
ind_sub <-
  purrr::map(summary$Samples, 1) %>%
  unlist() %>%
  as.character()

# Add ind_sub as column in summary to merge it later
summary$Sample <- ind_sub

# Make sure that all singletons from TreeCluster are kept
tokeep_singletons <-
  metadata %>%
  filter(`tc_0.003` == -1) %>%
  pull(Sample) %>%
  as.character()

# Make list of all individuals to keep in the tree
tokeep <- unique(c(tokeep_singletons, ind_sub))

# Make final list of ind to drop
todrop <-
  metadata %>%
  filter(!Sample %in% tokeep) %>%
  pull(Sample) %>%
  as.character()

# Make reduced metadata --------------------------------------------------------
metadata_red <-
  metadata %>%
  filter(!Sample %in% todrop) %>%
  left_join(., summary[,c("Sample", "short_haplogroup","n","n_ASW", "n_ESN", "n_ACB", "n_LWK", "n_MSL", "n_GWD", "n_YRI")])


# Match clusters across different groupings ------------------------------------

# Match/rename TreeCluster and RhierBAPS clusters based on rhb_01
# names(metadata_red)[6:12]
# "rhb_01"   "rhb_02"   "rhb_03"   "tc_0.006" "tc_0.005" "tc_0.004" "tc_0.003"

metadata_new <-
  metadata_red %>%
  filter(Sample != "RSRS") %>%
  match_clusters(., 6, 7, c(6:12)) %>%
  match_clusters(., 7, 8, c(6:12)) %>%
  match_clusters(., 8, 9, c(6:12)) %>%
  match_clusters(., 9, 10, c(6:12)) %>%
  match_clusters(., 10, 11, c(6:12)) %>%
  match_clusters(., 11, 12, c(6:12))

# Save table with renamed clusters for easier deciding on colors
metadata_new %>%
  write_csv(path_out_tree_summary_metadata)

################################################################################
# Match original and new cluster names
################################################################################

metadata_original <-
  metadata %>%
  rename_with(~ paste("original", .x, sep = "_")) %>%
  rename(
    Sample = original_Sample,
    Population = original_Population,
    Haplogroup = original_Haplogroup
  )

# Make summary table with original and matched names
metadata_new_with_original <- left_join(metadata_new, metadata_original)

# Left join column by column
all_metadata <-
  metadata_original %>%
  ljc(., metadata_new_with_original, "rhb_01") %>%
  ljc(., metadata_new_with_original, "rhb_02") %>%
  ljc(., metadata_new_with_original, "rhb_03") %>%
  ljc(., metadata_new_with_original, "tc_0.006") %>%
  ljc(., metadata_new_with_original, "tc_0.005") %>%
  ljc(., metadata_new_with_original, "tc_0.004") %>%
  ljc(., metadata_new_with_original, "tc_0.003") %>%
  ljc(., metadata_new_with_original, "SC") %>%
  ljc(., metadata_new_with_original, "SCL")

# Adjust RSRS
all_metadata[all_metadata$Sample == "RSRS", -1] <- NA
all_metadata[all_metadata$Sample == "RSRS", "Population"] <- "mt-MRCA"
all_metadata[all_metadata$Sample == "RSRS", "SC"] <- "mt-MRCA"
all_metadata[all_metadata$Sample == "RSRS", "SCL"] <- "mt-MRCA"
all_metadata[all_metadata$Sample == "RSRS", "Haplogroup"] <- "mt-MRCA"

# Save updated table with original and new cluster names
all_metadata %>% 
  write_csv(path_out_all_metadata)
