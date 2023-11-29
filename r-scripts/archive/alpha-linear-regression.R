## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(lme4)

load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

anyNA(tax_table(A_rumen_counts)[,"Phylum"])
anyNA(tax_table(H_rumen_counts)[,"Phylum"])

A_rumen_counts <- tax_fix(A_rumen_counts)
H_rumen_counts <- tax_fix(H_rumen_counts)

A_rumen_counts <- A_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_counts <- H_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

adivHr <- data.frame(
"Observed" = phyloseq::estimate_richness(H_rumen_counts, measures = "Observed"),
"Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
"Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment,
"Hour" = phyloseq::sample_data(H_rumen_counts)$Hour,
"Steer.ID" = phyloseq::sample_data(H_rumen_counts)$Steer.ID
)
head(adivHr)
write.csv(adivHr,"tables/holstein-rumen.csv", row.names = FALSE)

rm.shannon.hr <- lmer(Shannon ~ Treatment*Hour + (1|Steer.ID), data = adivHr)
summary(rm.shannon.hr)

glm.shannon.hr <- glm(Shannon ~ Treatment*Hour, data = adivHr)
summary(glm.shannon.hr)
