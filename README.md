[![DOI](https://zenodo.org/badge/679718318.svg)](https://zenodo.org/doi/10.5281/zenodo.10156922)

# mtDNA “Nomenclutter” and its Consequences on the Interpretation of Genetic Data

This repository includes scripts and data necessary to reproduce analyses and figures from the study "mtDNA “Nomenclutter” and its Consequences on the Interpretation of Genetic Data".

**Developers:** 
Vladimir Bajić and Vanessa Hava Schulmann

**Affiliation:**
[Human Biology and Primate Evolution, Freie Universität Berlin, Berlin, Germany](https://www.bcp.fu-berlin.de/en/biologie/arbeitsgruppen/zoologie/ag_nowick/index.html)

------------------------------------------------------------------------------------------------------------------------

# Data download

| Type      | File                               | url                                                                                                             | 
|-----------|------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| Metadata  | `igsr_sample.tsv`                  | https://www.internationalgenome.org/data-portal/sample                                                          | 
| Metadata  | `igsr_populations.tsv`             | https://www.internationalgenome.org/data-portal/population                                                      |
| Metadata  | `20140625_related_individuals.txt` | http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/20140625_related_individuals.txt                     |
| mtDNA Ref | `rCRS.fasta.gz`                    | http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/MT/rCRS.fasta.gz                          |
| mtDNA Ref | `RSRS.fasta.gz`                    | https://www.cell.com/cms/10.1016/j.ajhg.2012.03.002/attachment/34f040af-c97e-4bfa-854b-d2238be91d2e/mmc2.zip    |
| Data      | `mtDNA_1kgp.fasta.gz`              | http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/MT/chrMT_sequences_2534.20160505.fasta.gz |

------------------------------------------------------------------------------------------------------------------------

# Making input fasta file

We downloaded fasta file with full mitochondrial sequences from the 1000 Genomes Project. Then we added to it rCRS and RSRS sequences and named it `data/1KGP_rCRS_RSRS.fasta`.

> Note:<br/> 
    The revised Cambridge Reference Sequence (**rCRS**) was added as a reference for multiple sequence alignment.<br/> 
    The Reconstructed Sapiens Reference Sequence (**RSRS**) for rooting the phylogenetic tree. 

------------------------------------------------------------------------------------------------------------------------

# Removing relatives

Remove one of the individuals from the known related pairs based on the information from `20140625_related_individuals.txt`.

```bash
Rscript --vanilla scripts/remove_relatives.R
```

Output:
`out/fastas/1KGP_rCRS_RSRS_norelatives.fasta`

------------------------------------------------------------------------------------------------------------------------

# Manualy removing poly-c region and replacing N and M characters

Positions known to interfere with multiple-sequence alignment were removed or altered manually using `UGENE`. The spacer in the rCRS at position 3107 was changed from an "N" to an indel. The same was done for the spacers in the RSRS at positions 523-524. The positions in the poly-c region 303–315 and 16,183–16,194 were removed in all sequences. To avoid errors in downstream analysis character "M" (signaling a base ambiguity between "A" and "C") at position 12091 in individual HG01747 was manually changed to "C", a character present at this position in the majority of sequences in the dataset (exceptions NA18633 and NA19327 that have "A" at this position), as well as in the closest haplotype to HG01747.

Output: 
`out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced.fasta`

------------------------------------------------------------------------------------------------------------------------

# Multiple sequence alignment (MSA)

MSA was performed using MAFFT v7.520.

```bash
mafft out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced.fasta > out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_mafft.fasta
```

------------------------------------------------------------------------------------------------------------------------

# Manual post-alignment indel correction

Following good practice recommendations, manual post-alignment base correction (especialy around the original 3107 rCRS spacer) was performed using Unipro UGENE, ensuring that all indels align in the same pattern.

Output: 
`out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_mafft_indelcorr_postmafft.fa`

------------------------------------------------------------------------------------------------------------------------

# Creating African subset

```bash
Rscript --vanilla scripts/create_AFR_subset_FASTA.R
```

Output:
`out/fastas/1KGP_AFR_RSRS.fasta`
`out/fastas/1KGP_AFR.fasta`

------------------------------------------------------------------------------------------------------------------------

# Phylogeny construction

mtDNA phylogenetic trees were constructed by Maximum Likelihood (ML) with IQ-TREE. The ModelFinder Plus setting was used to determine the best-fitting model by calculating the Bayesian Information Criterion (BIC) and choosing the model which minimizes the BIC. Phylogenetic trees were rooted using the RSRS. 

```bash
iqtree -s out/fastas/1KGP_AFR_RSRS.fasta \
       -m MFP \
       -merit BIC \
       -o RSRS \
       -bb 1000 \
       -alrt 1000 \
       -nt AUTO \
       -prefix out/trees/1KGP_AFR_RSRS_iqtree
```

Output: `out/trees/1KGP_AFR_RSRS_iqtree.contree`

------------------------------------------------------------------------------------------------------------------------

# mtDNA groupings

## Haplogrep3
Haplogroup calling was performed using Haplogrep3 (v.3.2.1) command line tool, with FASTA as the input file, and the PhyloTree17 - Forensic Update (rCRS Human mtDNA) Version 1.2 (phylotree-fu-rcrs@1.2) as the tree. 

```bash
# create fasta files without gaps
tr -d '-' < out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced.fasta > out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_noGaps.fasta

# haplogrep3
haplogrep3 classify \
--in out/fastas/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_noGaps.fasta \
--tree phylotree-fu-rcrs@1.2 \
--out out/haplogrep3/1KGP_rCRS_RSRS_norelatives_nopolyc_NMreplaced_noGaps-haplogrep3_phylotree-fu-rcrs-1.2
```

## Nomenclature-Based Groupings (NBG)
Based on the output of HaploGrep3 we created two commonly occurring secondary NBGs: Single Character grouping (SC) and Single Character and L with one digit grouping (SCL) in R (see `scripts/make_all_metadata.R`).

## Algorithm-Based Groupings (ABG)
To obtain mitochondrial groupings based on mtDNA sequence similarity independent of traditional nomenclature, we performed two ABGs: rhierBAPS and TreeCluster.

### rhierBAPS (rhb)
The rhierBAPS groupings (rhb) were performed directly on MSA FASTA using the R package `rhierBAPS`. The rhierBAPS analysis was run with `keep.singletons=TRUE`, `max.depth=10`, and the default value for `n.pop`. For visualization purposes, we focussed on the first 3 levels of the rhierBAPS output.

```bash
Rscript --vanilla scripts/rhierBAPS.R
```

### TreeCluster (tc)
The TreeCluster groupings (tc) were performed on the mtDNA consensus trees outputted by the `IQ-TREE` using the command line program `TreeCluster` with eight different threshold values (`-t`) ranging from 0.001 to 0.008, and method (`-m`) set to its default: `max_clade`. The Max Clade method of TreeCluster clusters the leaves of the provided phylogenetic tree ensuring that the maximum pairwise distance between leaves in the cluster is at most equal to the specified threshold. For visualization purposes, we focussed on the output of TreeCluster with threshold values 0.003-0.006.

```bash
# TreeCluster
INPUT_FILE="out/trees/1KGP_AFR_RSRS_iqtree.contree"

for threshold in {0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008} 
    do TreeCluster.py -i $INPUT_FILE -t $threshold > out/treecluster/$(basename $INPUT_FILE)_t_$threshold.txt
done
```

------------------------------------------------------------------------------------------------------------------------

# Frequency-based analysis

## Frequency bar plots
Frequency bar plots based on different mitochondrial groupings were created in R using packages `tidyverse` and `ggpubr` (See `scripts/Fig_4_barplots.R`).

## Correspondence analysis (CA)
CA was performed using the R package `FactoMineR` and visualized using the R package `factoextra` (See `scripts/Fig_6_CA.R`).

------------------------------------------------------------------------------------------------------------------------

# Pairwise distance-based analyses
Pairwise distances between individuals from the same mitochondrial group were calculated in R using the packages `ape` and `tidyverse`.

## Multidimensional scaling (MDS)
MDS analysis was performed on the mtDNA pairwise distance matrix using the function `cmdscale()` from the R package `stats` and visualized using `ggplot2` (See `scripts/Fig_2_MDS_and_Map.R`).

## Violin plots
Violin plots of the mtDNA pairwise distances between individuals belonging to the same clusters produced by NBGs (SC and SCL), and ABGs (rhierBAPS and TreeCluster) were visualized using `ggplot2` and `ggpubr` (See `scripts/Fig_5_violinplots.R`).

------------------------------------------------------------------------------------------------------------------------

# Geographic map
The geographic map with populations was created using R packages `rworldmap` and `tidyverse` (See `scripts/Fig_2_MDS_and_Map.R`).

------------------------------------------------------------------------------------------------------------------------

# Phylogeny visualization
To obtain a schematic reduced phylogenetic tree (Figure 3) results from secondary nomenclature-based groupings (NBG) and algorithm-based groupings (ABG) were combined to identify a minimal set of representative samples that capture all unique combinations of different groupings across individuals. The phylogenetic tree was visualized using the R package `ggtree` (See `scripts/Fig_2_MDS_and_Map.R`).

------------------------------------------------------------------------------------------------------------------------

# Metadata
Collect all results and metadata into single file that will be used for ploting results.

```bash
# Metadata
Rscript --vanilla scripts/make_all_metadata.R
```
------------------------------------------------------------------------------------------------------------------------

# Plots

```bash
# Figure 2
Rscript --vanilla scripts/Fig_2_MDS_and_Map.R

# Figure 3
Rscript --vanilla scripts/Fig_3_summary_tree.R

# Figure 4
Rscript --vanilla scripts/Fig_4_barplots.R

# Figure 5
Rscript --vanilla scripts/Fig_5_violinplots.R

# Figure 6
Rscript --vanilla scripts/Fig_6_CA.R
```

------------------------------------------------------------------------------------------------------------------------
