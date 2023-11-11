################################################################################
# Script to create fasta file with AFR individuals
#
# This script:
#   - subsets 1KGP fasta file to individuals from AFR superpopulation
#   - creates a new fasta file for the subsetted data
#
# Author:
#   Vladimir Bajić
#
################################################################################


# Packages ---------------------------------------------------------------------
library(tidyverse)
library(seqinr)

# Set variables ----------------------------------------------------------------
path_metadata <- "data/igsr_samples.tsv"
path_fasta <- "out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_mafft_indelcorr_postmafft.fasta"
path_out_afr_fasta <- "out/fastas/1KGP_AFR.fasta"
path_out_afr_fasta_rsrs <- "out/fastas/1KGP_AFR_RSRS.fasta"

# Load metadata ----------------------------------------------------------------
# Subset to AFR Population
# downloaded from: https://www.internationalgenome.org/data-portal/sample
afr_samples <-
  read_tsv(path_metadata) %>%
  rename(Sample =  "Sample name", Superpopulation = "Superpopulation code") %>%
  select(Sample, Superpopulation)  %>%
  filter(Superpopulation == "AFR")  %>%
  pull(Sample)

# Load fasta -------------------------------------------------------------------
fasta <- read.fasta(path_fasta)

# Subset fasta based on sample IDs ---------------------------------------------
fasta_afr_subset_rsrs <- fasta[names(fasta) %in% c("RSRS", afr_samples)]
fasta_afr_subset <- fasta_afr_subset_rsrs[names(fasta_afr_subset_rsrs) != "RSRS"]

# Create fasta subset without RSRS ---------------------------------------------
write.fasta(sequences = fasta_afr_subset,
            names = names(fasta_afr_subset),
            file.out = path_out_afr_fasta)

# Create fasta subset with RSRS ------------------------------------------------
write.fasta(sequences = fasta_afr_subset_rsrs,
            names = names(fasta_afr_subset_rsrs),
            file.out = path_out_afr_fasta_rsrs)
