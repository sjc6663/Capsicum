# Check for Breed Difference, PERMANOVA
# 11-14-2022 SAB

## ---- setup ----
library(phyloseq)
library(microViz)

## ---- load relative abundance phyloseqs ----
load("ps-obj/phyloseq-fecal-samples-only-relabund.RData")
load("ps-obj/phyloseq-rumen-samples-only-relabund.RData")

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

## ---- PERMANOVA ----

# transform
transform_fecal <- microbiome::transform(fecal_rel, 'clr')
transform_rumen <- microbiome::transform(rumen_rel, 'clr')

# generate distance matrix
dist_matrix_fec <- phyloseq::distance(transform_fecal, method = "euclidean")
dist_matrix_rum <- phyloseq::distance(transform_rumen, method = "euclidean")

#ADONIS test
vegan::adonis2(dist_matrix_fec ~ phyloseq::sample_data(transform_fecal)$Breed)
vegan::adonis2(dist_matrix_rum ~ phyloseq::sample_data(transform_rumen)$Breed)
