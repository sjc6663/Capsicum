# Repeated Measures - Alpha Diversity Rumen Time Analysis
# April 11, 2023 - SAB

# setup ----
# load phyloseq objects
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)

# create data frame ----
adivH <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment,
  "Hour" = phyloseq::sample_data(H_rumen_counts)$Hour,
  "SteerID" = phyloseq::sample_data(H_rumen_counts)$Steer.ID
)

# lineplot ----

adivH$"Hour" <- factor(adivH$"Hour", 
                       levels = c("H0", "H2", "H6", "H12", "H18"))

adivH$"Treatment" <- factor(adivH$"Treatment", 
                            levels = c("Control", "RPC5", "RPC10", "RPC15"))

library(ggpubr)

line <- adivH %>% 
  ggline(x = "Hour", y = "Shannon", group = "Treatment", color = "Treatment", palette = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1"), add = "mean_sd", error.plot = "errorbar")

my_comparisons <- list( c("H0, H2"), c("H0, H6"), c("H0", "H12"), c("H0", "H18"), c("H2, H6"), c("H2", "H12"), c("H2", "H18"), c("H6", "H12"), c("H6", "H18"), c("H12", "H18"))

line + stat_compare_means(method = "anova", label.y = 1.4) + 
  stat_compare_means(aes(group = Treatment), label = "p.signif")

ggsave(filename = "plots/geomline-plot-Hrumen-time.jpeg", dpi = 600)
