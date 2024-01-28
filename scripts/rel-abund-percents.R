## Relative Abundance Percentages
# SAB 11-29-2023

## ---- setup ----
# load packages
library(BiocManager)
library(phyloseq)
library(vegan)
library(patchwork)
library(pairwiseAdonis)
library(microViz)
library(ggplot2)

# load phyloseq objects
load("ps-obj/phyloseq-fecal-samples-holstein-relabund.RData")

# check if we have NAs
anyNA(tax_table(H_fecal_rel)[,"Phylum"])

# tax fix our phyloseq object
H_fecal_rel <- tax_fix(H_fecal_rel)

H_fecal_rel <- H_fecal_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(H_fecal_rel, "Treatment")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Genus") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- dplyr::select(phy, c("Genus", "Sample", "Abundance"))

phy

### Select Rows we Want
grep("Limosilactobacillus", phy$Genus)
grep("Clostridium sensu stricto 18", phy$Genus)
phy2 <- phy[c(5, 10, 19, 36, 396, 397, 688, 689), ]

phy2
phy2 <- setorder(phy2, Genus)

phy2$percent <- phy2$Abundance * 100

phy2
phy2$percent <- round(phy2$percent, digits = 0)

phy2
