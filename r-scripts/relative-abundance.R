## Taxa Relative Abundance
# SAB 11/9/2022

## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(plyr)

# load phyloseq objects
psrel <- readRDS("ps-decontam-relabund.rds")
pscount <- readRDS("ps-decontam-filtered-counts.rds")

psrel

# merge everything at the Genus level so we don't have duplicates
psgrel <- tax_glom(psrel, "Genus")
psgcount <- tax_glom(pscount, "Genus")

## ---- separate out Rumen and Fecal Samples from relative abundance set ----
fecal_rel <- subset_samples(
  psgrel, 
  Sample.Type == "Fecal"
)

fecal_rel
sample_data(fecal_rel)

rumen_rel <- subset_samples(
  psgrel,
  Sample.Type == "Rumen"
)

rumen_rel
sample_data(rumen_rel)

# save these separate phyloseq objects
save(fecal_rel, file = "beta-diversity/phyloseq-fecal-samples-only-relabund.RData")
save(rumen_rel, file = "beta-diversity/phyloseq-rumen-samples-only-relabund.RData")

## ---- Tax Fix ----
### From here on out you have to run doubles, so one for fecal and one for rumen ###


# check if we have NAs
anyNA(tax_table(fecal_rel)[,"Phylum"])
anyNA(tax_table(rumen_rel)[,"Phylum"])

# tax fix our phyloseq object
fecal_rel <- tax_fix(fecal_rel)
rumen_rel <- tax_fix(rumen_rel)

fecal_rel <- fecal_rel %>% tax_fix(unknowns = c("Incertae Sedis"))
rumen_rel <- rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

## ---- Plot in microViz ----
# Fecal
fp_rel <- comp_barplot(fecal_rel,
                       tax_level = "Phylum",
                       group_by = "Treatment",
                       facet_by = "Breed")
fp_rel

# Rumen
rp_rel <- comp_barplot(rumen_rel,
                       tax_level = "Phylum",
                       group_by = "Treatment",
                       facet_by = "Breed")

rp_rel

