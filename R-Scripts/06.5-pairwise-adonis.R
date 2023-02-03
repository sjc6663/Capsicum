# Pairwise Adonis Test
# Holstein Fecal Samples
# SAB 01-05-2023

# setup
library(phyloseq)
library(vegan)
library(microViz)
library(BiocManager)

# install pairwise adonis
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

# load phyloseq objects
load("beta-diversity/phyloseq-fecal-samples-holstein-relabund.RData")

anyNA(tax_table(H_fecal_rel)[,"Phylum"])

H_fecal_rel <- tax_fix(H_fecal_rel)

H_fecal_rel <- H_fecal_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

# transform 
H_transform_fecal <- microbiome::transform(H_fecal_rel, 'clr')

# generate distance matrix
hol_fecal_dist_matrix <- phyloseq::distance(H_fecal_rel, method = "euclidean")

# pairwise adonis test
pairwise.adonis(hol_fecal_dist_matrix, sample_data(H_fecal_rel)$Treatment)
