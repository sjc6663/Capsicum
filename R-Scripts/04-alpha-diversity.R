# Alpha Diversity - Capsicum (Shannon's Diversity)
# SAB 11-02-2022

## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggpubr)
library(patchwork)

## ---- load phyloseq objects ----
load("alpha-diversity/phyloseq-fecal-samples-angus-counts.RData")
load("alpha-diversity/phyloseq-fecal-samples-holstein-counts.RData")
load("alpha-diversity/phyloseq-rumen-samples-angus-counts.RData")
load("alpha-diversity/phyloseq-rumen-samples-holstein-counts.RData")

# check if we have any NAs
anyNA(tax_table(A_rumen_counts)[,"Phylum"])
anyNA(tax_table(H_rumen_counts)[,"Phylum"])
anyNA(tax_table(A_fecal_counts)[,"Phylum"])
anyNA(tax_table(H_fecal_counts)[,"Phylum"])

# tax fix our phyloseq object
A_fecal_counts <- tax_fix(A_fecal_counts)
H_fecal_counts <- tax_fix(H_fecal_counts)
A_rumen_counts <- tax_fix(A_rumen_counts)
H_rumen_counts <- tax_fix(H_rumen_counts)

# individual tax fix
A_fecal_counts <- A_fecal_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_fecal_counts <- H_fecal_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
A_rumen_counts <- A_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_counts <- H_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

## ---- create data frames for each breed and type ----
adivHf <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_fecal_counts)$Treatment 
)
head(adivHf)
write.csv(adivHf,"alpha-diversity/tables/holstein-fecal.csv", row.names = FALSE)


adivHr <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment 
)
head(adivHr)
write.csv(adivHr,"alpha-diversity/tables/holstein-rumen.csv", row.names = FALSE)

adivAf <- data.frame(
  # "Observed" = phyloseq::estimate_richness(A_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_fecal_counts)$Treatment 
)
head(adivAf)
write.csv(adivAf,"alpha-diversity/tables/angus-fecal.csv", row.names = FALSE)

adivAr <- data.frame(
  # "Observed" = phyloseq::estimate_richness(A_rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_rumen_counts)$Treatment 
)
head(adivAr)
write.csv(adivAr,"alpha-diversity/tables/angus-rumen.csv", row.names = FALSE)

## ---- Shannon ANOVA with Holstein and Angus ----
testHf <- aov(Shannon ~ Treatment, data = adivHf)
summary(testHf)

testHr <- aov(Shannon ~ Treatment, data = adivHr)
summary(testHr)

testAf <- aov(Shannon ~ Treatment, data = adivAf)
summary(testAf)

testAr <- aov(Shannon ~ Treatment, data = adivAr)
summary(testAr)

## ---- BoxPlots ----
sample_data(H_fecal_counts)$"Treatment" <- factor(sample_data(H_fecal_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))
A <- plot_richness(H_fecal_counts, x="Treatment", measures=c("Shannon"), title = "A", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  theme_classic() +
  theme(legend.position = "none")
A
ggsave(plot = A, filename = "Plots/alpha-diversity-holstein-fecal.pdf", dpi = 600)

sample_data(H_rumen_counts)$"Treatment" <- factor(sample_data(H_rumen_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))
B <- plot_richness(H_rumen_counts, x="Treatment", measures=c("Shannon"), title = "B", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) + 
  theme_classic() +
  theme(legend.position = "none")
ggsave(plot = B, filename = "Plots/alpha-diversity-holstein-rumen.pdf", dpi = 600)

sample_data(A_fecal_counts)$"Treatment" <- factor(sample_data(A_fecal_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))
C <- plot_richness(A_fecal_counts, x="Treatment", measures=c("Shannon"), title = "C", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) + 
  theme_classic() +
  theme(legend.position = "none")
ggsave(plot = C, filename = "Plots/alpha-diversity-angus-fecal.pdf", dpi = 600)

sample_data(A_rumen_counts)$"Treatment" <- factor(sample_data(A_rumen_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))
D <- plot_richness(A_rumen_counts, x="Treatment", measures=c("Shannon"), title = "D", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) + 
  theme_classic() +
  theme(legend.position = "none")
ggsave(plot = D, filename = "Plots/alpha-diversity-angus-rumen.pdf", dpi = 600)

(A|B)/(C|D)
ggsave(filename = "Plots/shannon-diversity-all.pdf", dpi = 600)
