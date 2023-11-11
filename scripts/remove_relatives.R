################################################################################
# Script to create fasta file without known relatives
#
# This script:
#   - subsets 1KGP fasta file to individuals without known relatives
#
# Author:
#   Vladimir Bajić
#
################################################################################


# Packages ---------------------------------------------------------------------
library(tidyverse)
library(seqinr)

# Set variables ----------------------------------------------------------------
path_relatives <- "data/20140625_related_individuals.txt"
path_fasta <- "data/1KGP_rCRS_RSRS.fasta"
path_out <- "out/fastas/1KGP_rCRS_RSRS_norelatives.fasta"

# Load relatives info ----------------------------------------------------------
relatives <- read_delim(path_relatives)

# Load fasta -------------------------------------------------------------------
fasta <- read.fasta(path_fasta)

# Subset fasta based on sample IDs ---------------------------------------------
fasta_no_relatives <- fasta[!(names(fasta) %in% relatives$Sample)]

# Create fasta subset without relatives ----------------------------------------
write.fasta(sequences = fasta_no_relatives,
            names = names(fasta_no_relatives),
            file.out = path_out)