################################################################################
# Script to run rhierBAPS for the mtDNA from AFR superpop from 1KGP
#
# This script:
#   - runs rhierBAPS
#
# Author:
#   Vladimir Bajić
#
################################################################################

# Packages ---------------------------------------------------------------------
library(tidyverse)
library(rhierbaps)

# Set variables ----------------------------------------------------------------
#path_afr_fasta <- "data/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced.fa"
path_afr_fasta <- "out/fastas/1KGP_AFR.fasta"
out_name <- "out/rhierBAPS/rhb_1KGP_AFR_level10_keep_singletons_T.csv"
out_name_log <- "out/rhierBAPS/rhb_1KGP_AFR_level10_keep_singletons_T.lml_log"
max_depth <- 10
set.seed(1234)

# Create necessary dirs where the output will be stored ------------------------
dir.create("out/rhierBAPS", recursive = TRUE, showWarnings = FALSE)

# Load fasta as SNP matrix -----------------------------------------------------
snp_matrix <- load_fasta(path_afr_fasta, keep.singletons = TRUE)

# Clustering with hierBAPS -----------------------------------------------------
hb_results <- hierBAPS(snp_matrix,
                       max.depth = max_depth,
                       quiet = FALSE,
                       n.cores = 10)

# Save results -----------------------------------------------------------------
hb_results$partition.df  %>%
    write_csv(file = out_name, , quote = "all")

# Save lml logs ----------------------------------------------------------------
save_lml_logs(hb_results, out_name_log)